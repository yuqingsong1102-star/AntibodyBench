"""antibench/cli.py — 统一命令行入口。

用法（任意阶段均可单独运行）：

  # 阶段 1：数据准备（调用现有 scripts/data_prep/ 脚本）
  python -m antibench data_prep

  # 阶段 2：运行模型（单个模型或全部）
  python -m antibench run --model boltzgen
  python -m antibench run --model all

  # 阶段 2.5：收集候选（标准化）
  python -m antibench collect

  # 阶段 3a：生成 benchmark_1_runtime（pre-eval）
  python -m antibench benchmark --part runtime

  # 阶段 3b：生成 benchmark_2_common_k（eval 输入）
  python -m antibench benchmark --part b2 [--k 10]

  # 阶段 3c：生成 benchmark_1_hits（post-eval，需先跑 evaluation）
  python -m antibench benchmark --part hits

  # 阶段 4：evaluation pipeline
  python -m antibench evaluate [--candidates-csv path] [--steps 1 2 3 4]

  # 阶段 5：analysis（生成图表）
  python -m antibench analyze
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from antibench.utils import PROJECT_ROOT, PathManager, get_logger, load_config

logger = get_logger("cli")


# ─────────────────────────────────────────────────────────────────────────────
# 各阶段入口函数
# ─────────────────────────────────────────────────────────────────────────────

def cmd_data_prep(args: argparse.Namespace) -> int:
    """调用现有 scripts/data_prep/prepare_native_inputs.py。"""
    script = PROJECT_ROOT / "scripts" / "data_prep" / "prepare_native_inputs.py"
    if not script.exists():
        logger.error(f"data_prep 脚本不存在: {script}")
        return 1
    extra = args.extra or []
    cmd = [sys.executable, str(script)] + extra
    logger.info(f"[data_prep] 运行: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    return result.returncode


def cmd_run(args: argparse.Namespace) -> int:
    """运行指定模型（或全部模型）的 Runner。"""
    from antibench.dataset import load_targets
    from antibench.runner import get_runner, get_all_runners

    pm = PathManager()
    dataset_cfg = load_config("dataset")
    benchmark_cfg = load_config("benchmark")

    targets = load_targets(dataset_config=dataset_cfg, path_manager=pm)
    if args.target_ids:
        filter_set = set(args.target_ids)
        targets = [t for t in targets if t.target_id in filter_set]
    logger.info(f"[run] 共 {len(targets)} 个靶点")

    model_arg: str = args.model
    if model_arg == "all":
        model_names = benchmark_cfg.get("models", ["boltzgen", "mber_open", "germinal", "rfantibody"])
        runners = [get_runner(m, pm) for m in model_names]
    else:
        runners = [get_runner(model_arg, pm)]

    all_ok = True
    for runner in runners:
        inputs = runner._build_inputs(targets)
        metas = runner.run_all(
            inputs,
            force=args.force,
            continue_on_error=not args.fail_fast,
        )
        failed = [m for m in metas if m.status == "failed"]
        if failed:
            logger.warning(
                f"[{runner.model_name}] {len(failed)}/{len(metas)} 个靶点失败"
            )
            all_ok = False

    return 0 if all_ok else 1


def cmd_collect(args: argparse.Namespace) -> int:
    """收集所有模型的候选，写入 candidates_manifest.csv。"""
    from antibench.dataset import load_targets
    from antibench.collection import collect_all

    pm = PathManager()
    dataset_cfg = load_config("dataset")
    benchmark_cfg = load_config("benchmark")

    targets = load_targets(dataset_config=dataset_cfg, path_manager=pm)
    model_names: list[str] = benchmark_cfg.get(
        "models", ["boltzgen", "mber_open", "germinal", "rfantibody"]
    )

    out_path = Path(args.out) if args.out else pm.candidates_manifest
    collect_all(
        targets=targets,
        model_names=model_names,
        path_manager=pm,
        out_path=out_path,
        include_failed=args.include_failed,
    )
    return 0


def cmd_benchmark(args: argparse.Namespace) -> int:
    """生成 benchmark 报告。"""
    from antibench.benchmark import (
        benchmark_1_runtime,
        benchmark_1_hits,
        benchmark_2_common_k,
    )
    from antibench.candidate import load_candidates

    pm = PathManager()
    benchmark_cfg = load_config("benchmark")
    model_names: list[str] = benchmark_cfg.get(
        "models", ["boltzgen", "mber_open", "germinal", "rfantibody"]
    )

    part: str = args.part

    if part == "runtime":
        candidates = load_candidates(pm.candidates_manifest)
        benchmark_1_runtime(
            model_names=model_names,
            candidates=candidates,
            path_manager=pm,
        )

    elif part == "b2":
        k = args.k or benchmark_cfg.get("k_common", 10)
        candidates_csv = Path(args.candidates_csv) if args.candidates_csv else None
        benchmark_2_common_k(
            candidates_csv=candidates_csv,
            k=int(k),
            path_manager=pm,
            benchmark_config=benchmark_cfg,
        )

    elif part == "hits":
        eval_csv = Path(args.eval_csv) if args.eval_csv else pm.eval_final_csv
        if not eval_csv.exists():
            logger.error(f"eval_final.csv 不存在: {eval_csv}，请先运行 evaluate 阶段")
            return 1
        benchmark_1_hits(
            eval_final_csv=eval_csv,
            model_names=model_names,
            path_manager=pm,
            benchmark_config=benchmark_cfg,
        )

    else:
        logger.error(f"未知 --part 值: {part}，支持 runtime / b2 / hits")
        return 1

    return 0


def _build_eval_cmd(
    python_prefix: list[str],
    candidates_csv: Path,
    out_dir: Path,
    device: str,
    steps: list[int],
    use_pyrosetta: bool,
) -> list[str]:
    """构造 evaluation.run_evaluation 的调用命令。
    python_prefix 例如 ['python'] 或 ['conda', 'run', '-n', 'chai', 'python']。
    使用 -m 模块调用，避免 sys.path 问题。
    """
    cmd = python_prefix + [
        "-m", "evaluation.run_evaluation",
        "--candidates-csv", str(candidates_csv),
        "--out-dir", str(out_dir),
        "--device", device,
    ]
    for s in steps:
        cmd += ["--step", str(s)]
    if use_pyrosetta:
        cmd.append("--use-pyrosetta")
    return cmd


def cmd_evaluate(args: argparse.Namespace) -> int:
    """运行 evaluation pipeline（适配新接口）。

    多环境调度策略（来自 configs/evaluation.yaml）：
    - step1_python: Step 1 (Chai-1) 使用的 python 命令，例如 'conda run -n chai python'
    - default_python: Step 2-4 使用的 python 命令，例如 '/home/user/env/bin/python'
    若同时请求 step1 和其他步骤，会分两次 subprocess 调用。
    """
    pm = PathManager()
    eval_cfg = load_config("evaluation")

    # 候选 CSV：优先用命令行参数，其次用 config，最后用默认 benchmark_2
    if args.candidates_csv:
        candidates_csv = Path(args.candidates_csv)
    elif eval_cfg.get("candidates_csv"):
        candidates_csv = PROJECT_ROOT / eval_cfg["candidates_csv"]
    else:
        k = load_config("benchmark").get("k_common", 10)
        candidates_csv = pm.benchmark_2_csv(int(k))

    if not candidates_csv.exists():
        logger.error(
            f"候选 CSV 不存在: {candidates_csv}\n"
            f"请先运行：python -m antibench benchmark --part b2"
        )
        return 1

    out_dir = Path(args.out_dir) if args.out_dir else pm.evaluation_dir
    steps: list[int] = args.steps or eval_cfg.get("steps", [1, 2, 3, 4])
    device: str = args.device or eval_cfg.get("device", "cuda:0")
    use_pyrosetta: bool = args.use_pyrosetta

    # ── 多环境调度 ────────────────────────────────────────────────────────
    # step1_python 示例："conda run -n chai python" 或 "/abs/path/python"
    step1_python_str: str = eval_cfg.get("step1_python", sys.executable)
    default_python_str: str = eval_cfg.get("default_python", sys.executable)
    step1_prefix = step1_python_str.split()      # 字符串拆为列表
    default_prefix = default_python_str.split()

    logger.info(f"[evaluate] 输入: {candidates_csv}")
    logger.info(f"[evaluate] 输出: {out_dir}")
    logger.info(f"[evaluate] 步骤: {steps}")
    logger.info(f"[evaluate] Step1 Python: {step1_python_str}")
    logger.info(f"[evaluate] 其他 Python: {default_python_str}")

    step1_steps = [s for s in steps if s == 1]
    other_steps = [s for s in steps if s != 1]

    # Step 1 单独用 chai 环境调用
    if step1_steps:
        cmd = _build_eval_cmd(step1_prefix, candidates_csv, out_dir, device, step1_steps, use_pyrosetta)
        logger.info(f"[evaluate/step1] 运行: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=str(PROJECT_ROOT))
        if result.returncode != 0:
            logger.error(f"[evaluate/step1] 失败 (exit={result.returncode})")
            if not other_steps:
                return result.returncode
            logger.warning("[evaluate/step1] 继续执行后续步骤...")

    # Step 2-4 用主环境调用
    if other_steps:
        cmd = _build_eval_cmd(default_prefix, candidates_csv, out_dir, device, other_steps, use_pyrosetta)
        logger.info(f"[evaluate/steps{other_steps}] 运行: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=str(PROJECT_ROOT))
        if result.returncode != 0:
            logger.error(f"[evaluate/steps{other_steps}] 失败 (exit={result.returncode})")
            return result.returncode

    return 0


def cmd_analyze(args: argparse.Namespace) -> int:
    """运行 analysis（生成图表）。

    代理给 analysis/analyze_results.py，补充合理默认值：
    - eval_dirs 默认 outputs/evaluation
    - out_dir   默认 outputs/analysis
    """
    pm = PathManager()
    script = PROJECT_ROOT / "analysis" / "analyze_results.py"
    if not script.exists():
        logger.error(f"分析脚本不存在: {script}")
        return 1

    # eval_dirs：可多值，默认 outputs/evaluation
    eval_dirs: list[str] = args.eval_dirs or [str(pm.evaluation_dir)]
    out_dir: str = args.out_dir or str(pm.analysis_dir)

    cmd = [
        sys.executable,
        str(script),
        "--eval-dirs", *eval_dirs,
        "--out-dir", out_dir,
    ]
    if args.model_labels:
        cmd += ["--model-labels", *args.model_labels]
    if args.top_k is not None:
        cmd += ["--top-k", str(args.top_k)]
    if args.rank_by:
        cmd += ["--rank-by", args.rank_by]

    logger.info(f"[analyze] eval_dirs={eval_dirs}")
    logger.info(f"[analyze] out_dir={out_dir}")
    logger.info(f"[analyze] 运行: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(PROJECT_ROOT))
    return result.returncode


# ─────────────────────────────────────────────────────────────────────────────
# 参数解析
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m antibench",
        description="AntibodyBench MVP pipeline 统一入口",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="stage", required=True)

    # ── data_prep ────────────────────────────────────────────────────────
    p = sub.add_parser("data_prep", help="数据准备（调用现有 data_prep 脚本）")
    p.add_argument("extra", nargs=argparse.REMAINDER, help="透传给 prepare_native_inputs.py 的额外参数")

    # ── run ──────────────────────────────────────────────────────────────
    p = sub.add_parser("run", help="运行模型")
    p.add_argument(
        "--model", required=True,
        help="模型名（boltzgen / mber_open / germinal / rfantibody / all）",
    )
    p.add_argument("--target-ids", nargs="*", help="只运行指定靶点（空格分隔）")
    p.add_argument("--force", action="store_true", help="忽略 .done 标记，强制重新运行")
    p.add_argument("--fail-fast", action="store_true", help="遇到第一个失败即中止")

    # ── collect ──────────────────────────────────────────────────────────
    p = sub.add_parser("collect", help="收集并标准化候选")
    p.add_argument("--out", help="输出 CSV 路径（默认：outputs/candidates/candidates_manifest.csv）")
    p.add_argument("--include-failed", action="store_true", help="包含状态为 failed 的候选")

    # ── benchmark ────────────────────────────────────────────────────────
    p = sub.add_parser("benchmark", help="生成 benchmark 报告")
    p.add_argument(
        "--part", required=True,
        choices=["runtime", "b2", "hits"],
        help="runtime=benchmark_1_runtime | b2=benchmark_2_common_k | hits=benchmark_1_hits",
    )
    p.add_argument("--k", type=int, help="benchmark_2 的 K 值（默认从 config 读取）")
    p.add_argument("--candidates-csv", help="候选 CSV 路径（b2 用）")
    p.add_argument("--eval-csv", help="eval_final.csv 路径（hits 用）")

    # ── evaluate ─────────────────────────────────────────────────────────
    p = sub.add_parser("evaluate", help="运行 evaluation pipeline")
    p.add_argument("--candidates-csv", help="候选 CSV 输入路径（默认：benchmark_2 输出）")
    p.add_argument("--out-dir", help="评估输出目录（默认：outputs/evaluation/）")
    p.add_argument("--steps", nargs="*", type=int, help="运行哪些步骤（1~4，默认全部）")
    p.add_argument("--device", help="计算设备（默认：cuda:0）")
    p.add_argument("--use-pyrosetta", action="store_true", help="启用 PyRosetta ddG")

    # ── analyze ──────────────────────────────────────────────────────────
    p = sub.add_parser(
        "analyze",
        help="生成图表与统计表格（代理 analysis/analyze_results.py）",
    )
    p.add_argument(
        "--eval-dirs", nargs="+", default=None,
        metavar="DIR",
        help="evaluation 结果目录（可多个，默认：outputs/evaluation）",
    )
    p.add_argument(
        "--model-labels", nargs="+", default=None,
        metavar="LABEL",
        help="与 --eval-dirs 对应的模型标签",
    )
    p.add_argument(
        "--out-dir", default=None,
        help="图表输出目录（默认：outputs/analysis）",
    )
    p.add_argument(
        "--top-k", type=int, default=None,
        help="每个靶点取 top-K 候选参与分析（默认 0=全部）",
    )
    p.add_argument(
        "--rank-by", default=None,
        choices=["chai1_iptm", "pdockq2", "composite"],
        help="排名指标（默认：chai1_iptm）",
    )

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    dispatch = {
        "data_prep": cmd_data_prep,
        "run": cmd_run,
        "collect": cmd_collect,
        "benchmark": cmd_benchmark,
        "evaluate": cmd_evaluate,
        "analyze": cmd_analyze,
    }
    return dispatch[args.stage](args)


if __name__ == "__main__":
    sys.exit(main())
