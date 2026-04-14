"""antibench/utils.py — 配置读取、路径管理、日志、计时。

所有模块从这里取配置和路径，不在各自脚本里硬编码。
"""
from __future__ import annotations

import csv
import json
import logging
import os
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any

# ── 可选依赖：yaml ────────────────────────────────────────────────────────
try:
    import yaml as _yaml
    _HAS_YAML = True
except ImportError:  # pragma: no cover
    _HAS_YAML = False


# ─────────────────────────────────────────────────────────────────────────────
# 项目根目录（utils.py 在 antibench/ 下，所以 parent.parent = 项目根）
# ─────────────────────────────────────────────────────────────────────────────
PROJECT_ROOT: Path = Path(__file__).parent.parent.resolve()


# ─────────────────────────────────────────────────────────────────────────────
# 配置加载
# ─────────────────────────────────────────────────────────────────────────────

def load_yaml(path: Path | str) -> dict[str, Any]:
    """读取 YAML 文件，返回 dict。不依赖 yaml 时抛出友好错误。"""
    if not _HAS_YAML:
        raise ImportError("请先安装 pyyaml：conda install pyyaml 或 pip install pyyaml")
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        return _yaml.safe_load(f) or {}


def load_config(name: str) -> dict[str, Any]:
    """加载配置文件，支持三种写法：

    * load_config('benchmark')              → configs/benchmark.yaml
    * load_config('benchmark.yaml')         → configs/benchmark.yaml
    * load_config('configs/benchmark.yaml') → configs/benchmark.yaml（项目根相对路径）
    * load_config('/abs/path/to/file.yaml') → 直接用绝对路径
    """
    p = Path(name)

    # 绝对路径：直接使用
    if p.is_absolute():
        if not p.exists():
            raise FileNotFoundError(f"配置文件不存在: {p}")
        return load_yaml(p)

    # 去掉 .yaml 后缀，再去掉 configs/ 前缀，得到裸名
    stem = name
    if stem.endswith(".yaml"):
        stem = stem[:-5]
    if stem.startswith("configs/"):
        stem = stem[len("configs/"):]
    # 支持 configs\\xxx（Windows 路径）
    if stem.startswith("configs\\"):
        stem = stem[len("configs\\"):]

    config_path = PROJECT_ROOT / "configs" / f"{stem}.yaml"
    if not config_path.exists():
        raise FileNotFoundError(
            f"配置文件不存在: {config_path}\n"
            f"  支持的写法：load_config('benchmark') / load_config('benchmark.yaml') / "
            f"load_config('configs/benchmark.yaml')"
        )
    return load_yaml(config_path)


def load_model_config(model_name: str) -> dict[str, Any]:
    """加载 configs/models/<model_name>.yaml。"""
    config_path = PROJECT_ROOT / "configs" / "models" / f"{model_name}.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"模型配置文件不存在: {config_path}")
    return load_yaml(config_path)


# ─────────────────────────────────────────────────────────────────────────────
# 路径管理
# ─────────────────────────────────────────────────────────────────────────────

class PathManager:
    """统一管理项目内各阶段的输出路径，避免路径字符串散落各处。"""

    def __init__(self, root: Path | str | None = None) -> None:
        self.root = Path(root).resolve() if root else PROJECT_ROOT
        self.outputs = self.root / "outputs"

    # ── 数据层 ─────────────────────────────────────────────────────────────
    @property
    def data_dir(self) -> Path:
        return self.root / "data"

    @property
    def native_inputs_dir(self) -> Path:
        return self.root / "native_inputs"

    # ── 输出层 ─────────────────────────────────────────────────────────────
    @property
    def data_prep_dir(self) -> Path:
        return self.outputs / "data_prep"

    @property
    def cropped_dir(self) -> Path:
        return self.data_prep_dir / "cropped"

    def run_dir(self, model_name: str, target_id: str | None = None) -> Path:
        base = self.outputs / "runs" / model_name
        return base / target_id if target_id else base

    @property
    def candidates_dir(self) -> Path:
        return self.outputs / "candidates"

    @property
    def candidates_manifest(self) -> Path:
        return self.candidates_dir / "candidates_manifest.csv"

    @property
    def benchmark_dir(self) -> Path:
        return self.outputs / "benchmark"

    def benchmark_2_csv(self, k: int) -> Path:
        return self.benchmark_dir / f"benchmark_2_k{k}.csv"

    @property
    def benchmark_1_runtime_csv(self) -> Path:
        return self.benchmark_dir / "benchmark_1_runtime.csv"

    @property
    def benchmark_1_hits_csv(self) -> Path:
        return self.benchmark_dir / "benchmark_1_hits.csv"

    @property
    def evaluation_dir(self) -> Path:
        return self.outputs / "evaluation"

    @property
    def eval_final_csv(self) -> Path:
        return self.evaluation_dir / "eval_final.csv"

    @property
    def analysis_dir(self) -> Path:
        return self.outputs / "analysis"

    def ensure(self, *paths: Path) -> None:
        """确保给定目录都存在。"""
        for p in paths:
            p.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# 日志
# ─────────────────────────────────────────────────────────────────────────────

def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """返回带格式的 logger，同一 name 不重复添加 handler。"""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter(
                fmt="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )
        logger.addHandler(handler)
    logger.setLevel(level)
    return logger


def setup_file_logger(logger: logging.Logger, log_path: Path) -> None:
    """给已有 logger 追加文件输出。"""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    logger.addHandler(fh)


# ─────────────────────────────────────────────────────────────────────────────
# 计时
# ─────────────────────────────────────────────────────────────────────────────

@contextmanager
def timer(label: str = "", logger: logging.Logger | None = None):
    """上下文管理器，测量代码块执行时间并打印/记录。"""
    t0 = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - t0
        msg = f"{label} 耗时 {elapsed:.2f}s" if label else f"耗时 {elapsed:.2f}s"
        if logger:
            logger.info(msg)
        else:
            print(msg)


# ─────────────────────────────────────────────────────────────────────────────
# 通用 IO
# ─────────────────────────────────────────────────────────────────────────────

def read_csv(path: Path) -> list[dict[str, str]]:
    """读取 CSV，返回 list[dict]。文件不存在时返回空列表。"""
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as f:
        return [{str(k): str(v) for k, v in row.items()} for row in csv.DictReader(f)]


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    """写入 CSV，自动创建父目录。fieldnames 为 None 时从第一行推断。"""
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fieldnames or list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> dict:
    """读取 JSON，文件不存在或解析失败时返回空 dict。"""
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def write_json(path: Path, data: dict) -> None:
    """写入 JSON，自动创建父目录。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def safe_str(value: Any) -> str:
    """转成字符串，None → 空字符串。"""
    return "" if value is None else str(value).strip()


def to_project_rel(path: Path | str | None, root: Path | None = None) -> str:
    """把绝对路径转成相对于项目根的字符串，用于 CSV 中存储可移植路径。"""
    if path is None:
        return ""
    root = root or PROJECT_ROOT
    try:
        return str(Path(path).resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)
