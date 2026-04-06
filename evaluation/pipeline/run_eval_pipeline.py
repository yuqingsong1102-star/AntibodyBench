#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from evaluation.ingest.build_external_eval_inputs import build_external_eval_inputs
from evaluation.judge.run_primary_judge import run_primary_judge
from evaluation.judge.run_secondary_judge import run_secondary_judge
from evaluation.metrics.compute_affinity_proxy import compute_affinity_proxy
from evaluation.metrics.compute_binding_plausibility import compute_binding_plausibility
from evaluation.metrics.compute_developability import compute_developability
from evaluation.metrics.compute_robustness import compute_robustness
from evaluation.report.build_model_comparison import build_model_comparison
from evaluation.report.build_scorecards import build_scorecards
from evaluation.report.visualize_external_metrics import generate_external_figures
from evaluation.report.write_external_report import write_external_report


def run_eval_pipeline(
  *,
  native_roots: list[Path],
  index_csv: Path,
  out_dir: Path,
  judge_manifest: Path | None,
  secondary_judge_manifest: Path | None,
  fallback_to_source_structure: bool,
  judge_name: str,
  n_requested_seeds: int,
  max_candidates_per_model_sample: int,
  judge_cmd: str | None,
) -> None:
  out_dir.mkdir(parents=True, exist_ok=True)
  candidates_csv = out_dir / "external_eval_candidates.csv"
  candidate_selection_audit_csv = out_dir / "external_eval_candidates_selection_audit.csv"
  seed_level_csv = out_dir / "external_eval_seed_level.csv"
  sequence_csv = out_dir / "external_eval_sequence_level.csv"
  model_summary_csv = out_dir / "external_eval_model_summary.csv"
  sample_summary_csv = out_dir / "external_eval_sample_summary.csv"
  fig_dir = out_dir / "figures"
  report_md = out_dir / "report.md"
  judge_input_dir = out_dir / "judge_inputs"
  primary_judge_run_dir = out_dir / "primary_judge_runs"

  build_external_eval_inputs(
    native_roots=native_roots,
    index_csv=index_csv,
    out_csv=candidates_csv,
    job_dir=judge_input_dir,
    judge_name=judge_name,
    n_requested_seeds=n_requested_seeds,
    max_candidates_per_model_sample=max_candidates_per_model_sample,
    selection_audit_csv=candidate_selection_audit_csv,
  )
  run_primary_judge(
    candidates_csv=candidates_csv,
    out_csv=seed_level_csv,
    judge_manifest=judge_manifest,
    fallback_to_source_structure=fallback_to_source_structure,
    judge_cmd=judge_cmd,
    run_dir=primary_judge_run_dir,
  )
  run_secondary_judge(seed_level_csv=seed_level_csv, out_csv=seed_level_csv, judge_manifest=secondary_judge_manifest)
  compute_binding_plausibility(in_csv=seed_level_csv, out_csv=seed_level_csv)
  compute_affinity_proxy(in_csv=seed_level_csv, out_csv=seed_level_csv)
  compute_developability(in_csv=seed_level_csv, out_csv=seed_level_csv)
  compute_robustness(in_csv=seed_level_csv, out_csv=seed_level_csv)
  build_scorecards(seed_level_csv=seed_level_csv, sequence_csv=sequence_csv)
  build_model_comparison(
    sequence_csv=sequence_csv,
    seed_level_csv=seed_level_csv,
    model_summary_csv=model_summary_csv,
    sample_summary_csv=sample_summary_csv,
  )
  generate_external_figures(sequence_csv=sequence_csv, model_summary_csv=model_summary_csv, fig_dir=fig_dir)
  write_external_report(
    model_summary_csv=model_summary_csv,
    seed_level_csv=seed_level_csv,
    report_path=report_md,
    fig_dir=fig_dir,
  )


def main() -> int:
  parser = argparse.ArgumentParser(description="运行统一外部评估流水线（按 UNIFIED_EXTERNAL_EVALUATION_PLAN）")
  parser.add_argument(
    "--native-root",
    action="append",
    type=Path,
    help="native 输出目录（可重复指定）",
  )
  parser.add_argument(
    "--index-csv",
    type=Path,
    default=Path("data/prepared/dataset_index_ready.csv"),
    help="样本索引 CSV",
  )
  parser.add_argument(
    "--out-dir",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark"),
    help="输出目录",
  )
  parser.add_argument(
    "--judge-manifest",
    type=Path,
    default=None,
    help="primary judge 结果 manifest，可选",
  )
  parser.add_argument(
    "--secondary-judge-manifest",
    type=Path,
    default=None,
    help="secondary judge 结果 manifest，可选",
  )
  parser.add_argument(
    "--fallback-to-source-structure",
    action="store_true",
    help="没有 judge manifest 时，回退使用模型原始结构作为占位输入",
  )
  parser.add_argument("--judge-name", type=str, default="chai1_primary", help="主 judge 名称")
  parser.add_argument("--n-requested-seeds", type=int, default=3, help="每条候选请求 seeds 数")
  parser.add_argument(
    "--max-candidates-per-model-sample",
    type=int,
    default=10,
    help="每个 model/sample 保留的 top-K 唯一候选序列数",
  )
  parser.add_argument(
    "--judge-cmd",
    type=str,
    default=None,
    help="直接执行主 judge 的命令模板，可用占位符如 __JOB_JSON__/__JUDGE_SEED__/__OUTPUT_DIR__/__RESULT_JSON__",
  )
  args = parser.parse_args()

  if args.native_root:
    native_roots = args.native_root
  elif Path("outputs/native_predictions_run2").exists():
    native_roots = [Path("outputs/native_predictions_run2")]
  else:
    native_roots = [Path("outputs/native_predictions_real"), Path("outputs/native_predictions_smoke")]

  run_eval_pipeline(
    native_roots=native_roots,
    index_csv=args.index_csv,
    out_dir=args.out_dir,
    judge_manifest=args.judge_manifest,
    secondary_judge_manifest=args.secondary_judge_manifest,
    fallback_to_source_structure=args.fallback_to_source_structure,
    judge_name=args.judge_name,
    n_requested_seeds=max(1, args.n_requested_seeds),
    max_candidates_per_model_sample=max(1, args.max_candidates_per_model_sample),
    judge_cmd=args.judge_cmd,
  )
  print("[OK] 外部评估流水线完成")
  print(f"[OK] 结果目录: {args.out_dir}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
