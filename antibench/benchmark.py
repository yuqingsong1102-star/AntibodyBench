"""antibench/benchmark.py — 两个 benchmark 层次。

benchmark_1_runtime（pre-eval）
  来源：RunMetadata（outputs/runs/<model>/<target_id>/run_meta.json）
  内容：GPU 时间统计、成功率、候选总数
  时机：模型运行完成后即可生成，无需 eval_final.csv

benchmark_1_hits（post-eval）
  来源：eval_final.csv ⊕ RunMetadata
  内容：命中率、time-to-first-hit、GPU hours per hit
  时机：evaluation pipeline 完成后生成

benchmark_2_common_k（eval 输入）
  来源：candidates_manifest.csv
  操作：每个 (model, target_id) 按 rank_in_model 取前 K_common 个
  输出：outputs/benchmark/benchmark_2_k{K}.csv（evaluation pipeline 的输入）
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from antibench.candidate import Candidate, RunMetadata, load_candidates
from antibench.utils import (
    PathManager,
    get_logger,
    load_config,
    read_csv,
    write_csv,
)

logger = get_logger("benchmark")


# ─────────────────────────────────────────────────────────────────────────────
# Benchmark 1 — runtime summary（pre-eval）
# ─────────────────────────────────────────────────────────────────────────────

RUNTIME_FIELDS = [
    "model_name",
    "target_id",
    "status",
    "wall_clock_sec",
    "gpu_count",
    "gpu_hours",
    "total_candidates",
    "exit_code",
    "error_summary",
    "run_dir",
]


def collect_run_metas(
    model_names: list[str],
    target_ids: list[str] | None = None,
    path_manager: PathManager | None = None,
) -> list[RunMetadata]:
    """从 outputs/runs/<model>/<target_id>/run_meta.json 加载所有 RunMetadata。"""
    pm = path_manager or PathManager()
    metas: list[RunMetadata] = []

    for model_name in model_names:
        model_run_dir = pm.run_dir(model_name)
        if not model_run_dir.exists():
            continue
        for target_dir in sorted(model_run_dir.iterdir()):
            if not target_dir.is_dir():
                continue
            target_id = target_dir.name
            if target_ids and target_id not in target_ids:
                continue
            meta_path = target_dir / "run_meta.json"
            if not meta_path.exists():
                continue
            try:
                import json
                data = json.loads(meta_path.read_text(encoding="utf-8"))
                metas.append(RunMetadata.from_dict(data))
            except Exception as exc:
                logger.warning(f"读取 run_meta.json 失败 {meta_path}: {exc}")

    logger.info(f"[benchmark] 加载了 {len(metas)} 条 RunMetadata")
    return metas


def _count_candidates_for(
    model_name: str,
    target_id: str,
    candidates: list[Candidate],
) -> int:
    return sum(
        1 for c in candidates
        if c.model_name == model_name and c.target_id == target_id
    )


def benchmark_1_runtime(
    model_names: list[str],
    candidates: list[Candidate] | None = None,
    path_manager: PathManager | None = None,
    out_path: Path | None = None,
) -> list[dict[str, str]]:
    """生成 benchmark_1_runtime.csv（pre-eval，不需要 eval 结果）。

    Args:
        model_names: 要统计的模型列表。
        candidates: 已收集的 Candidate 列表，用于统计候选数。可为 None。
        path_manager: PathManager 实例。
        out_path: 输出路径，None 则用默认。

    Returns:
        list of row dicts（同时写入 CSV）。
    """
    pm = path_manager or PathManager()
    out_path = out_path or pm.benchmark_1_runtime_csv

    metas = collect_run_metas(model_names, path_manager=pm)
    candidates = candidates or []

    rows: list[dict[str, str]] = []
    for meta in metas:
        n_candidates = _count_candidates_for(
            meta.model_name, meta.target_id, candidates
        )
        row = {
            "model_name": meta.model_name,
            "target_id": meta.target_id,
            "status": meta.status,
            "wall_clock_sec": f"{meta.wall_clock_sec:.2f}" if meta.wall_clock_sec is not None else "",
            "gpu_count": str(meta.gpu_count),
            "gpu_hours": f"{meta.gpu_hours:.4f}" if meta.gpu_hours is not None else "",
            "total_candidates": str(n_candidates),
            "exit_code": str(meta.exit_code) if meta.exit_code is not None else "",
            "error_summary": meta.error_summary,
            "run_dir": str(meta.run_dir or ""),
        }
        rows.append(row)

    write_csv(out_path, rows, RUNTIME_FIELDS)
    logger.info(f"[benchmark_1_runtime] 已写入 {len(rows)} 行到 {out_path}")

    # 打印摘要
    _print_runtime_summary(rows, model_names)
    return rows


def _print_runtime_summary(rows: list[dict[str, str]], model_names: list[str]) -> None:
    # 使用 rows 中实际存在的 model_name，而非外部传入的 model_names 列表
    # （外部列表可能来自 config，名称可能与实际目录扫描结果不一致）
    actual_models = list(dict.fromkeys(r["model_name"] for r in rows if r.get("model_name")))
    if not actual_models:
        logger.warning("[benchmark_1_runtime] 未找到任何运行元数据，请检查 outputs/runs/ 目录")
        return
    print("\n─── benchmark_1_runtime 摘要 ───")
    for model in actual_models:
        model_rows = [r for r in rows if r["model_name"] == model]
        ok = sum(1 for r in model_rows if r["status"] == "ok")
        total = len(model_rows)
        total_gpu_h = sum(
            float(r["gpu_hours"]) for r in model_rows if r.get("gpu_hours")
        )
        print(
            f"  {model:<15} 完成率={ok}/{total}  "
            f"总GPU时间={total_gpu_h:.2f}h"
        )
    print()


# ─────────────────────────────────────────────────────────────────────────────
# Benchmark 1 — hits summary（post-eval）
# ─────────────────────────────────────────────────────────────────────────────

HITS_FIELDS = [
    "model_name",
    "target_id",
    "k_evaluated",
    "hit_count",
    "hit_rate",
    "first_hit_rank",
    "wall_clock_sec",
    "gpu_hours",
    "gpu_hours_per_hit",
]


def _is_hit(row: dict[str, str], thresholds: dict[str, float]) -> bool:
    """按 configs/benchmark.yaml 中的 hit_thresholds 判断一行是否为 hit。
    全部满足才算 hit（AND 逻辑）。
    """
    for field_name, threshold in thresholds.items():
        value_str = row.get(field_name, "")
        if not value_str:
            return False
        try:
            if float(value_str) < threshold:
                return False
        except ValueError:
            return False
    return True


def benchmark_1_hits(
    eval_final_csv: Path,
    model_names: list[str] | None = None,
    path_manager: PathManager | None = None,
    benchmark_config: dict | None = None,
    out_path: Path | None = None,
) -> list[dict[str, str]]:
    """生成 benchmark_1_hits.csv（post-eval，需要 eval_final.csv）。

    Args:
        eval_final_csv: evaluation pipeline 产出的 eval_final.csv 路径。
        model_names: 要统计的模型列表，None 则自动从 CSV 中读取。
        path_manager: PathManager 实例。
        benchmark_config: configs/benchmark.yaml，None 时自动加载。
        out_path: 输出路径。

    Returns:
        list of row dicts。
    """
    pm = path_manager or PathManager()
    out_path = out_path or pm.benchmark_1_hits_csv
    bcfg = benchmark_config or load_config("benchmark")
    hit_thresholds: dict[str, float] = {
        k: float(v) for k, v in bcfg.get("hit_thresholds", {}).items()
    }

    eval_rows = read_csv(eval_final_csv)
    if not eval_rows:
        logger.warning(f"eval_final.csv 为空或不存在: {eval_final_csv}")
        return []

    # 检查 hit_thresholds 所需字段是否存在
    if hit_thresholds and eval_rows:
        for field in hit_thresholds:
            non_empty = [r for r in eval_rows if r.get(field, "").strip()]
            if not non_empty:
                logger.warning(
                    f"[benchmark_1_hits] 字段 '{field}' 在 eval_final.csv 中全部为空！\n"
                    f"  → hit 将全部为 0。请确保 Step 1 (Chai-1) 已成功运行，\n"
                    f"    或检查 evaluation/steps/chai1_judge.py 是否正常产出该字段。"
                )

    # 加载 RunMetadata 用于 GPU 时间
    all_model_names = model_names or sorted({r.get("model_name", "") for r in eval_rows} - {""})
    metas_by_key: dict[tuple[str, str], RunMetadata] = {}
    for meta in collect_run_metas(all_model_names, path_manager=pm):
        metas_by_key[(meta.model_name, meta.target_id)] = meta

    # 按 (model, target_id) 分组
    from collections import defaultdict
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in eval_rows:
        model = row.get("model_name", "")
        target = row.get("target_id", "")
        if model and target:
            groups[(model, target)].append(row)

    output_rows: list[dict[str, str]] = []
    for (model, target_id), group_rows in sorted(groups.items()):
        if model_names and model not in model_names:
            continue

        # 按 rank_in_model 排序（数字小 = 内部排名高）
        try:
            group_rows = sorted(group_rows, key=lambda r: int(r.get("rank_in_model") or 9999))
        except ValueError:
            pass

        hits = [r for r in group_rows if _is_hit(r, hit_thresholds)]
        hit_count = len(hits)
        hit_rate = hit_count / len(group_rows) if group_rows else 0.0
        first_hit_rank = int(hits[0].get("rank_in_model") or 0) if hits else -1

        meta = metas_by_key.get((model, target_id))
        gpu_hours = meta.gpu_hours if meta else None
        gpu_hours_per_hit = (
            gpu_hours / hit_count if gpu_hours is not None and hit_count > 0 else None
        )

        output_rows.append({
            "model_name": model,
            "target_id": target_id,
            "k_evaluated": str(len(group_rows)),
            "hit_count": str(hit_count),
            "hit_rate": f"{hit_rate:.4f}",
            "first_hit_rank": str(first_hit_rank),
            "wall_clock_sec": f"{meta.wall_clock_sec:.2f}" if meta and meta.wall_clock_sec is not None else "",
            "gpu_hours": f"{gpu_hours:.4f}" if gpu_hours is not None else "",
            "gpu_hours_per_hit": f"{gpu_hours_per_hit:.4f}" if gpu_hours_per_hit is not None else "",
        })

    write_csv(out_path, output_rows, HITS_FIELDS)
    logger.info(f"[benchmark_1_hits] 已写入 {len(output_rows)} 行到 {out_path}")
    return output_rows


# ─────────────────────────────────────────────────────────────────────────────
# Benchmark 2 — common K 截取
# ─────────────────────────────────────────────────────────────────────────────

def benchmark_2_common_k(
    candidates: list[Candidate] | None = None,
    candidates_csv: Path | None = None,
    k: int | None = None,
    path_manager: PathManager | None = None,
    benchmark_config: dict | None = None,
    out_path: Path | None = None,
) -> tuple[list[Candidate], Path]:
    """按统一 K 截取候选，输出 benchmark_2_k{K}.csv 作为 evaluation 输入。

    每个 (model, target_id) 取 rank_in_model ≤ K 的候选。
    rank_in_model 由 collection.py 按模型内部得分排定，只用于本模型内部截取。

    Args:
        candidates: Candidate 列表，None 时从 CSV 加载。
        candidates_csv: 候选 CSV 路径，candidates 为 None 时读取。
        k: 每 (model, target_id) 取多少个候选，None 则从 config 读取。
        path_manager: PathManager 实例。
        benchmark_config: configs/benchmark.yaml，None 时自动加载。
        out_path: 输出 CSV 路径。

    Returns:
        (截取后的 Candidate 列表, 输出 CSV 路径)
    """
    pm = path_manager or PathManager()
    bcfg = benchmark_config or load_config("benchmark")
    k = k if k is not None else int(bcfg.get("k_common", 10))
    out_path = out_path or pm.benchmark_2_csv(k)

    if candidates is None:
        csv_path = candidates_csv or pm.candidates_manifest
        candidates = load_candidates(csv_path)

    if not candidates:
        logger.warning("[benchmark_2] 没有候选可截取")
        return [], out_path

    # 按 (model, target_id) 分组，取 rank_in_model <= k（已按 rank_in_model 升序）
    from collections import defaultdict
    groups: dict[tuple[str, str], list[Candidate]] = defaultdict(list)
    for c in candidates:
        if not c.is_reference:  # 天然参考不参与 K 截取
            groups[(c.model_name, c.target_id)].append(c)

    selected: list[Candidate] = []
    for (model, target_id), group in groups.items():
        group_sorted = sorted(group, key=lambda c: c.rank_in_model)
        selected.extend(group_sorted[:k])

    # 追加所有 native_reference（不受 K 限制）
    references = [c for c in candidates if c.is_reference]
    selected.extend(references)

    logger.info(
        f"[benchmark_2] K={k}，选出 {len(selected)} 个候选 "
        f"（其中 {len(references)} 个 native_reference）"
    )

    from antibench.candidate import save_candidates
    save_candidates(selected, out_path)
    logger.info(f"[benchmark_2] 已写入到 {out_path}")

    # 打印各模型摘要
    _print_b2_summary(selected, k)
    return selected, out_path


def _print_b2_summary(candidates: list[Candidate], k: int) -> None:
    from collections import defaultdict
    counts: dict[str, int] = defaultdict(int)
    for c in candidates:
        if not c.is_reference:
            counts[c.model_name] += 1

    print(f"\n─── benchmark_2_common_k (K={k}) 摘要 ───")
    for model, count in sorted(counts.items()):
        print(f"  {model:<15} {count} 个候选进入评估")
    print()
