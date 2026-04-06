#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from evaluation.report.build_model_comparison import build_model_comparison
from evaluation.report.build_scorecards import build_scorecards


def build_external_summaries(
  *,
  seed_csv: Path,
  sequence_csv: Path,
  model_summary_csv: Path,
  sample_summary_csv: Path,
) -> tuple[Path, Path, Path]:
  out_seq = build_scorecards(seed_level_csv=seed_csv, sequence_csv=sequence_csv)
  out_model, out_sample = build_model_comparison(
    sequence_csv=out_seq,
    seed_level_csv=seed_csv,
    model_summary_csv=model_summary_csv,
    sample_summary_csv=sample_summary_csv,
  )
  return out_seq, out_model, out_sample


def main() -> int:
  parser = argparse.ArgumentParser(description="兼容入口：从 seed-level 主表构建汇总")
  parser.add_argument(
    "--seed-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="seed 级主表",
  )
  parser.add_argument(
    "--sequence-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_sequence_level.csv"),
    help="sequence-level 汇总表",
  )
  parser.add_argument(
    "--model-summary-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_model_summary.csv"),
    help="模型汇总表",
  )
  parser.add_argument(
    "--sample-summary-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_sample_summary.csv"),
    help="样本汇总表",
  )
  args = parser.parse_args()
  sequence_csv, model_csv, sample_csv = build_external_summaries(
    seed_csv=args.seed_csv,
    sequence_csv=args.sequence_csv,
    model_summary_csv=args.model_summary_csv,
    sample_summary_csv=args.sample_summary_csv,
  )
  print(f"[OK] candidate 汇总: {sequence_csv}")
  print(f"[OK] 模型汇总: {model_csv}")
  print(f"[OK] 样本汇总: {sample_csv}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())