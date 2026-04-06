#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from evaluation.metrics.compute_affinity_proxy import compute_affinity_proxy
from evaluation.metrics.compute_binding_plausibility import compute_binding_plausibility
from evaluation.metrics.compute_developability import compute_developability
from evaluation.metrics.compute_robustness import compute_robustness


def compute_phase1_scorecard(in_csv: Path, out_csv: Path) -> Path:
  # 兼容入口：顺序执行新分层评分模块，最终仍输出 seed-level 表。
  compute_binding_plausibility(in_csv=in_csv, out_csv=out_csv)
  compute_affinity_proxy(in_csv=out_csv, out_csv=out_csv)
  compute_developability(in_csv=out_csv, out_csv=out_csv)
  compute_robustness(in_csv=out_csv, out_csv=out_csv)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="兼容入口：运行 Phase1 seed-level 评分卡")
  parser.add_argument(
    "--in-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="seed 级主表",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="输出 seed 级主表（可覆盖输入）",
  )
  args = parser.parse_args()
  out_csv = compute_phase1_scorecard(in_csv=args.in_csv, out_csv=args.out_csv)
  print(f"[OK] seed-level 评分卡: {out_csv}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())