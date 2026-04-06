#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from evaluation.ingest.build_external_eval_inputs import build_external_eval_inputs


def build_external_eval_candidates(
  *,
  native_roots: list[Path],
  index_csv: Path,
  out_csv: Path,
  job_dir: Path,
  max_candidates_per_model_sample: int = 10,
  selection_audit_csv: Path | None = None,
) -> Path:
  return build_external_eval_inputs(
    native_roots=native_roots,
    index_csv=index_csv,
    out_csv=out_csv,
    job_dir=job_dir,
    max_candidates_per_model_sample=max_candidates_per_model_sample,
    selection_audit_csv=selection_audit_csv,
  )


def main() -> int:
  parser = argparse.ArgumentParser(description="兼容入口：构建外部评估候选表与 judge 输入模板")
  parser.add_argument("--native-root", action="append", type=Path, help="native 输出目录（可重复指定）")
  parser.add_argument("--index-csv", type=Path, default=Path("data/prepared/dataset_index_ready.csv"), help="样本索引 CSV")
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_candidates.csv"),
    help="候选表输出路径",
  )
  parser.add_argument(
    "--job-dir",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/judge_inputs"),
    help="judge 输入模板目录",
  )
  parser.add_argument(
    "--max-candidates-per-model-sample",
    type=int,
    default=10,
    help="每个 model/sample 保留的 top-K 唯一候选序列数",
  )
  parser.add_argument(
    "--selection-audit-csv",
    type=Path,
    default=None,
    help="候选筛选审计表输出路径，默认与候选表同目录",
  )
  args = parser.parse_args()

  native_roots = args.native_root or [Path("outputs/native_predictions_run2")]
  out_csv = build_external_eval_candidates(
    native_roots=native_roots,
    index_csv=args.index_csv,
    out_csv=args.out_csv,
    job_dir=args.job_dir,
    max_candidates_per_model_sample=max(1, args.max_candidates_per_model_sample),
    selection_audit_csv=args.selection_audit_csv,
  )
  print(f"[OK] 候选表: {out_csv}")
  print(f"[OK] 候选筛选审计表: {args.selection_audit_csv or args.out_csv.with_name(f'{args.out_csv.stem}_selection_audit.csv')}")
  print(f"[OK] judge 输入目录: {args.job_dir}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())