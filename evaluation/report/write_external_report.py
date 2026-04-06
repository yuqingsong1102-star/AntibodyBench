#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _load_csv(path: Path) -> list[dict[str, str]]:
  if not path.exists():
    return []
  with path.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def write_external_report(
  *,
  model_summary_csv: Path,
  seed_level_csv: Path,
  report_path: Path,
  fig_dir: Path,
) -> Path:
  model_rows = _load_csv(model_summary_csv)
  seed_rows = _load_csv(seed_level_csv)
  report_path.parent.mkdir(parents=True, exist_ok=True)

  judge_modes = sorted({(row.get("judge_mode") or "").strip() for row in seed_rows if (row.get("judge_mode") or "").strip()})
  lines: list[str] = []
  lines.append("# External Binding Evaluation Report (Phase 1)")
  lines.append("")
  lines.append("## Scope")
  lines.append("- Official focus: viability, binding plausibility, affinity proxy, robustness, and developability.")
  lines.append("- Native-similarity metrics are intentionally excluded from this report.")
  if judge_modes:
    lines.append(f"- Judge modes present: {', '.join(judge_modes)}")
  lines.append("")
  lines.append("## Scorecard Layers")
  lines.append("1. Viability")
  lines.append("2. Binding plausibility")
  lines.append("3. Affinity proxy")
  lines.append("4. Robustness")
  lines.append("5. Developability")
  lines.append("")
  lines.append("## Key Figures")
  for name in [
    "fig_01_viability_rate.png",
    "fig_02_bindability_pass_at_k.png",
    "fig_03_pdockq2_boxplot.png",
    "fig_04_external_iptm_boxplot.png",
    "fig_05_interface_bsa_dg_scatter.png",
    "fig_06_surface_hydrophobicity_binder_score_scatter.png",
    "fig_07_robust_pass_rate.png",
  ]:
    lines.append(f"- `{fig_dir / name}`")
  lines.append("")
  lines.append("## Summary By Model")
  lines.append("")
  lines.append("| model | n_candidates | viability_rate | pass@1 | pass@5 | pass@10 | strong_binder_proxy_rate | robust_pass_rate | median_pdockq2 | median_external_iptm | median_interface_bsa | median_interface_dG |")
  lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
  for row in model_rows:
    lines.append(
      "| {model} | {n_candidates} | {viability_rate} | {bindability_pass_at_1} | {bindability_pass_at_5} | {bindability_pass_at_10} | {strong_binder_proxy_rate} | {robust_pass_rate} | {median_pdockq2} | {median_external_iptm} | {median_interface_bsa} | {median_interface_dG} |".format(
        model=row.get("model", ""),
        n_candidates=row.get("n_candidates", ""),
        viability_rate=row.get("viability_rate", ""),
        bindability_pass_at_1=row.get("bindability_pass_at_1", ""),
        bindability_pass_at_5=row.get("bindability_pass_at_5", ""),
        bindability_pass_at_10=row.get("bindability_pass_at_10", ""),
        strong_binder_proxy_rate=row.get("strong_binder_proxy_rate", ""),
        robust_pass_rate=row.get("robust_pass_rate", ""),
        median_pdockq2=row.get("median_pdockq2", ""),
        median_external_iptm=row.get("median_external_iptm", ""),
        median_interface_bsa=row.get("median_interface_bsa", ""),
        median_interface_dG=row.get("median_interface_dG", ""),
      )
    )
  lines.append("")
  lines.append("## Notes")
  lines.append("- Phase 1 uses a unified external scorecard; imported/fallback judge outputs are both explicitly marked by `judge_mode`.")
  lines.append("- `pass@K` is defined on candidate-level binding pass, not on similarity to native antibody.")
  lines.append("- Under single-seed fallback, robustness is a provisional signal and should not be over-interpreted.")
  lines.append("- Report intentionally excludes TM-score / DockQ / CDR-H3 RMSD rankings.")

  report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
  return report_path


def main() -> int:
  parser = argparse.ArgumentParser(description="写出 Phase1 外部评估报告")
  parser.add_argument(
    "--model-summary-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_model_summary.csv"),
    help="模型汇总表",
  )
  parser.add_argument(
    "--seed-level-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="seed 级主表",
  )
  parser.add_argument(
    "--report-path",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/report.md"),
    help="报告输出路径",
  )
  parser.add_argument(
    "--fig-dir",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/figures"),
    help="图表目录",
  )
  args = parser.parse_args()
  out = write_external_report(
    model_summary_csv=args.model_summary_csv,
    seed_level_csv=args.seed_level_csv,
    report_path=args.report_path,
    fig_dir=args.fig_dir,
  )
  print(f"[OK] 报告: {out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())