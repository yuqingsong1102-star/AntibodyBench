"""antibench/runner.py — 统一 Runner 接口 + 4 个模型的薄包装。

设计原则：
- BaseRunner 定义统一接口，每个模型只需覆盖少量方法
- 实际执行委托给现有 scripts/native_runners/*.sh（不重写模型逻辑）
- RunMetadata 自动记录：开始/结束时间、exit code、日志路径、状态
- 支持断点续跑：检测 .done 文件跳过已完成靶点
- Python 是主流程编排层，shell 是最薄的执行层
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from antibench.candidate import RunMetadata
from antibench.dataset import Target
from antibench.utils import (
    PROJECT_ROOT,
    PathManager,
    get_logger,
    load_model_config,
    write_json,
)

logger = get_logger("runner")


# ─────────────────────────────────────────────────────────────────────────────
# ModelInput — runner 的输入对象
# ─────────────────────────────────────────────────────────────────────────────

class ModelInput:
    """单个靶点对单个模型的输入描述。"""
    def __init__(
        self,
        target: Target,
        model_name: str,
        input_dir: Path,
        model_config: dict[str, Any],
    ) -> None:
        self.target = target
        self.model_name = model_name
        self.input_dir = input_dir           # native_inputs/<model_subdir>/<target_id>/
        self.model_config = model_config


# ─────────────────────────────────────────────────────────────────────────────
# BaseRunner
# ─────────────────────────────────────────────────────────────────────────────

class BaseRunner(ABC):
    """所有模型 Runner 的基类。"""

    def __init__(
        self,
        model_config: dict[str, Any] | None = None,
        path_manager: PathManager | None = None,
    ) -> None:
        self.model_name: str = self._model_name()
        self.cfg = model_config or load_model_config(self.model_name)
        self.pm = path_manager or PathManager()

    @staticmethod
    @abstractmethod
    def _model_name() -> str:
        """返回模型名称字符串（与 configs/models/<name>.yaml 一致）。"""

    # ── 断点续跑 ─────────────────────────────────────────────────────────

    def _done_flag(self, target_id: str) -> Path:
        done_name = self.pm.root / "configs" / "benchmark.yaml"
        # 从 benchmark 配置读取 done 文件名
        done_filename = ".done"
        try:
            import yaml
            cfg = yaml.safe_load(done_name.read_text())
            done_filename = cfg.get("done_filename", ".done")
        except Exception:
            pass
        return self.pm.run_dir(self.model_name, target_id) / done_filename

    def is_done(self, target_id: str) -> bool:
        """检测该靶点是否已成功完成（存在 .done 文件）。"""
        return self._done_flag(target_id).exists()

    def _mark_done(self, target_id: str) -> None:
        flag = self._done_flag(target_id)
        flag.parent.mkdir(parents=True, exist_ok=True)
        flag.touch()

    # ── 元数据持久化 ─────────────────────────────────────────────────────

    def _meta_path(self, target_id: str) -> Path:
        return self.pm.run_dir(self.model_name, target_id) / "run_meta.json"

    def _write_meta(self, meta: RunMetadata) -> None:
        """将 RunMetadata 写入 run_meta.json。"""
        write_json(self._meta_path(meta.target_id), meta.to_dict())

    def _load_meta(self, target_id: str) -> RunMetadata | None:
        path = self._meta_path(target_id)
        if not path.exists():
            return None
        data = json.loads(path.read_text(encoding="utf-8"))
        return RunMetadata.from_dict(data)

    # ── 主运行接口 ───────────────────────────────────────────────────────

    def run(self, model_input: ModelInput, *, force: bool = False) -> RunMetadata:
        """运行单个靶点。

        Args:
            model_input: 包含 Target 和模型配置的输入对象。
            force: True 时忽略 .done 标记，强制重新运行。

        Returns:
            RunMetadata，无论成功与否都返回。
        """
        target_id = model_input.target.target_id

        if not force and self.is_done(target_id):
            logger.info(f"[{self.model_name}] 跳过（已完成）: {target_id}")
            existing = self._load_meta(target_id)
            return existing or RunMetadata(
                model_name=self.model_name,
                target_id=target_id,
                status="skipped",
            )

        run_dir = self.pm.run_dir(self.model_name, target_id)
        run_dir.mkdir(parents=True, exist_ok=True)
        log_path = run_dir / "run.log"

        meta = RunMetadata(
            model_name=self.model_name,
            target_id=target_id,
            start_time=datetime.now(tz=timezone.utc),
            gpu_count=int(self.cfg.get("gpu_count", 1)),
            log_path=log_path,
            run_dir=run_dir,
            status="running",
        )
        self._write_meta(meta)

        logger.info(f"[{self.model_name}] 开始运行: {target_id}")
        t0 = time.perf_counter()

        try:
            exit_code = self._execute(model_input, run_dir=run_dir, log_path=log_path)
        except Exception as exc:
            exit_code = -1
            meta.error_summary = str(exc)
            logger.error(f"[{self.model_name}] 运行异常 {target_id}: {exc}")

        elapsed = time.perf_counter() - t0
        meta.end_time = datetime.now(tz=timezone.utc)
        meta.wall_clock_sec = elapsed
        meta.exit_code = exit_code

        if exit_code == 0:
            meta.status = "ok"
            self._mark_done(target_id)
            logger.info(f"[{self.model_name}] 完成: {target_id} ({elapsed:.1f}s)")
        else:
            meta.status = "failed"
            logger.warning(f"[{self.model_name}] 失败 (exit={exit_code}): {target_id}")

        self._write_meta(meta)
        return meta

    def run_all(
        self,
        inputs: list[ModelInput],
        *,
        force: bool = False,
        continue_on_error: bool = True,
    ) -> list[RunMetadata]:
        """串行运行所有靶点，返回 RunMetadata 列表。"""
        results: list[RunMetadata] = []
        for i, inp in enumerate(inputs, 1):
            logger.info(
                f"[{self.model_name}] 进度 {i}/{len(inputs)}: {inp.target.target_id}"
            )
            meta = self.run(inp, force=force)
            results.append(meta)
            if meta.status == "failed" and not continue_on_error:
                logger.error(f"[{self.model_name}] 中止（continue_on_error=False）")
                break
        return results

    # ── 子类实现 ─────────────────────────────────────────────────────────

    @abstractmethod
    def _execute(
        self,
        model_input: ModelInput,
        *,
        run_dir: Path,
        log_path: Path,
    ) -> int:
        """执行实际运行逻辑，返回 exit code。子类必须实现。"""

    # ── 内部工具 ─────────────────────────────────────────────────────────

    def _run_shell(
        self,
        cmd: list[str],
        *,
        env: dict[str, str] | None = None,
        log_path: Path | None = None,
        timeout_sec: int | None = None,
    ) -> int:
        """执行 shell 命令，将 stdout/stderr 写入 log_path。"""
        merged_env = {**os.environ, **(env or {})}
        log_sink = open(log_path, "w", encoding="utf-8") if log_path else None
        try:
            proc = subprocess.run(
                cmd,
                stdout=log_sink,
                stderr=subprocess.STDOUT,
                env=merged_env,
                timeout=timeout_sec,
            )
            return proc.returncode
        except subprocess.TimeoutExpired:
            logger.warning(f"命令超时 ({timeout_sec}s): {' '.join(cmd)}")
            return 124  # 与 bash timeout 一致
        except OSError as exc:
            logger.error(f"命令执行失败: {exc}")
            return 1
        finally:
            if log_sink:
                log_sink.close()

    def _build_inputs(self, targets: list[Target]) -> list[ModelInput]:
        """为给定靶点列表构造 ModelInput 列表。"""
        input_subdir = self.cfg.get("input_subdir", self.model_name)
        inputs = []
        for target in targets:
            input_dir = self.pm.native_inputs_dir / input_subdir / target.target_id
            inputs.append(
                ModelInput(
                    target=target,
                    model_name=self.model_name,
                    input_dir=input_dir,
                    model_config=self.cfg,
                )
            )
        return inputs


# ─────────────────────────────────────────────────────────────────────────────
# 具体 Runner 实现（薄包装，委托给现有 shell）
# ─────────────────────────────────────────────────────────────────────────────

class _ShellRunner(BaseRunner):
    """通用 shell-based Runner：调用 scripts/native_runners/<model>.sh。"""

    def _execute(
        self,
        model_input: ModelInput,
        *,
        run_dir: Path,
        log_path: Path,
    ) -> int:
        runner_script = PROJECT_ROOT / self.cfg["runner_script"]
        if not runner_script.exists():
            logger.error(f"runner 脚本不存在: {runner_script}")
            return 1

        # 传入 runner shell 所需的标准参数
        # （与 scripts/run.sh 约定的接口一致）
        cmd = [
            "bash",
            str(runner_script),
            "--sample-id", model_input.target.target_id,
            "--sample-input-dir", str(model_input.input_dir),
            "--sample-output-dir", str(run_dir),
            "--project-root", str(PROJECT_ROOT),
        ]

        # 注入模型专属环境变量（runner shell 通过环境变量读取配置）
        env = self._build_env(model_input)

        return self._run_shell(cmd, env=env, log_path=log_path)

    def _build_env(self, model_input: ModelInput) -> dict[str, str]:
        """构造传给 runner shell 的环境变量。子类可 override 追加。"""
        env: dict[str, str] = {}
        workdir = self.cfg.get("workdir", "")
        if workdir:
            model_key = self.model_name.upper().replace("-", "_")
            env[f"{model_key}_WORKDIR"] = str(workdir)
        return env


class BoltzGenRunner(_ShellRunner):
    @staticmethod
    def _model_name() -> str:
        return "boltzgen"


class MBEROpenRunner(_ShellRunner):
    @staticmethod
    def _model_name() -> str:
        return "mber_open"


class GerminalRunner(_ShellRunner):
    @staticmethod
    def _model_name() -> str:
        return "germinal"


class RFAntibodyRunner(_ShellRunner):
    @staticmethod
    def _model_name() -> str:
        return "rfantibody"


# ─────────────────────────────────────────────────────────────────────────────
# Runner 工厂
# ─────────────────────────────────────────────────────────────────────────────

_RUNNER_REGISTRY: dict[str, type[BaseRunner]] = {
    "boltzgen": BoltzGenRunner,
    "mber_open": MBEROpenRunner,
    "mber-open": MBEROpenRunner,
    "germinal": GerminalRunner,
    "rfantibody": RFAntibodyRunner,
    "RFantibody": RFAntibodyRunner,
}


def get_runner(
    model_name: str,
    path_manager: PathManager | None = None,
) -> BaseRunner:
    """根据模型名称返回对应的 Runner 实例。"""
    cls = _RUNNER_REGISTRY.get(model_name)
    if cls is None:
        supported = sorted(set(_RUNNER_REGISTRY.keys()))
        raise ValueError(
            f"不支持的模型: {model_name!r}。支持的模型: {supported}"
        )
    return cls(path_manager=path_manager)


def get_all_runners(path_manager: PathManager | None = None) -> list[BaseRunner]:
    """返回所有模型的 Runner 实例（去重）。"""
    seen: set[type] = set()
    runners = []
    for cls in _RUNNER_REGISTRY.values():
        if cls not in seen:
            seen.add(cls)
            runners.append(cls(path_manager=path_manager))
    return runners
