#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from build_evaluation_long import build_evaluation_long
from compute_cdr_h3_rmsd import compute_cdr_h3_rmsd
from compute_sequence_metrics import compute_sequence_metrics
from visualize_metrics import main as visualize_main


def run_pipeline(
  *,
  native_roots: list[Path],
  index_csv: Path,
  out_dir: Path,
) -> None:
  out_dir.mkdir(parents=True, exist_ok=True)
  long_csv = out_dir / "evaluation_long.csv"
  long_seq_csv = out_dir / "evaluation_long_seq.csv"
  long_metrics_csv = out_dir / "evaluation_long_metrics.csv"
  summary_model_csv = out_dir / "summary_by_model.csv"
  summary_sample_csv = out_dir / "summary_by_sample.csv"
  fig_dir = out_dir / "figures"
  report_md = out_dir / "report.md"

  build_evaluation_long(native_roots=native_roots, index_csv=index_csv, out_csv=long_csv)
  compute_sequence_metrics(
    in_csv=long_csv,
    out_csv=long_seq_csv,
    summary_by_model_csv=summary_model_csv,
    summary_by_sample_csv=summary_sample_csv,
  )
  compute_cdr_h3_rmsd(in_csv=long_seq_csv, out_csv=long_metrics_csv)
  compute_sequence_metrics(
    in_csv=long_metrics_csv,
    out_csv=long_metrics_csv,
    summary_by_model_csv=summary_model_csv,
    summary_by_sample_csv=summary_sample_csv,
  )

  # 复用 visualize 脚本 CLI 逻辑，避免重复实现。
  import sys

  old_argv = sys.argv
  try:
    sys.argv = [
      "visualize_metrics.py",
      "--metrics-csv",
      str(long_metrics_csv),
      "--summary-csv",
      str(summary_model_csv),
      "--fig-dir",
      str(fig_dir),
      "--report-path",
      str(report_md),
    ]
    visualize_main()
  finally:
    sys.argv = old_argv


def main() -> int:
  parser = argparse.ArgumentParser(description="一键运行轨道B半自动评估流水线")
  parser.add_argument(
    "--native-root",
    action="append",
    type=Path,
    help="native 输出目录（可重复指定），默认 real + smoke",
  )
  parser.add_argument(
    "--index-csv",
    type=Path,
    default=Path("inputs/antibody_datasets/dataset_index_h3_annotated.csv"),
    help="索引 CSV",
  )
  parser.add_argument(
    "--out-dir",
    type=Path,
    default=Path("outputs/evaluation/all_models"),
    help="评估结果目录",
  )
  args = parser.parse_args()

  native_roots = args.native_root or [
    Path("outputs/native_predictions_real"),
    Path("outputs/native_predictions_smoke"),
  ]
  run_pipeline(native_roots=native_roots, index_csv=args.index_csv, out_dir=args.out_dir)

  print("[OK] 轨道B评估流水线完成")
  print(f"[OK] 结果目录: {args.out_dir}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

