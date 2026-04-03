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


def write_report(summary_csv: Path, report_path: Path, fig_dir: Path | None = None) -> Path:
  summary_rows = _load_csv(summary_csv)
  report_path.parent.mkdir(parents=True, exist_ok=True)

  figure_names = [
    "fig_01_success_rate.png",
    "fig_02_cdr_h3_rmsd_boxplot.png",
    "fig_03_runtime_boxplot.png",
    "fig_04_sample_model_heatmap.png",
    "fig_05_pareto_runtime_vs_rmsd.png",
  ]

  lines: list[str] = []
  lines.append("# Evaluation Report (All Models)")
  lines.append("")
  lines.append("## Key Figures")
  for name in figure_names:
    if fig_dir is None:
      lines.append(f"- `{name}`")
    else:
      lines.append(f"- `{fig_dir / name}`")
  lines.append("")
  lines.append("## Summary By Model")
  lines.append("")
  lines.append("| model | n_total | n_success | success_rate | median_runtime_sec | median_cdr_h3_rmsd | median_tm_score | median_dockq_score | median_irmsd |")
  lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
  for r in summary_rows:
    lines.append(
      "| {model} | {n_total} | {n_success} | {success_rate} | {median_runtime_sec} | {median_cdr_h3_rmsd} | {median_tm_score} | {median_dockq_score} | {median_irmsd} |".format(
        model=r.get("model", ""),
        n_total=r.get("n_total", ""),
        n_success=r.get("n_success", ""),
        success_rate=r.get("success_rate", ""),
        median_runtime_sec=r.get("median_runtime_sec", ""),
        median_cdr_h3_rmsd=r.get("median_cdr_h3_rmsd", ""),
        median_tm_score=r.get("median_tm_score", ""),
        median_dockq_score=r.get("median_dockq_score", ""),
        median_irmsd=r.get("median_irmsd", ""),
      )
    )
  lines.append("")
  lines.append("## Notes")
  lines.append("- 失败样本已保留在评估长表中，不会被静默过滤。")
  lines.append("- 当结构结果缺失时，结构相关指标会保留为空字符串。")

  report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
  return report_path


def main() -> int:
  parser = argparse.ArgumentParser(description="从模型汇总表生成 Markdown 报告")
  parser.add_argument(
    "--summary-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/summary_by_model.csv"),
    help="模型汇总表",
  )
  parser.add_argument(
    "--report-path",
    type=Path,
    default=Path("outputs/evaluation/all_models/report.md"),
    help="报告输出路径",
  )
  parser.add_argument(
    "--fig-dir",
    type=Path,
    default=Path("outputs/evaluation/all_models/figures"),
    help="图像目录（用于写入报告中的文件名）",
  )
  args = parser.parse_args()

  out = write_report(summary_csv=args.summary_csv, report_path=args.report_path, fig_dir=args.fig_dir)
  print(f"[OK] 报告: {out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())