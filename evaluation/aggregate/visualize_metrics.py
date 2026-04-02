#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import io
import math
import contextlib
from pathlib import Path

_MPL_ERR = ""
try:
  with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore
except Exception as e:  # pragma: no cover
  plt = None
  _MPL_ERR = f"{type(e).__name__}: {e}"


_MINIMAL_PNG_1X1 = (
  b"\x89PNG\r\n\x1a\n"
  b"\x00\x00\x00\rIHDR"
  b"\x00\x00\x00\x01\x00\x00\x00\x01\x08\x06\x00\x00\x00"
  b"\x1f\x15\xc4\x89"
  b"\x00\x00\x00\nIDATx\x9cc`\x00\x00\x00\x02\x00\x01"
  b"\xe2!\xbc3"
  b"\x00\x00\x00\x00IEND\xaeB`\x82"
)


def _write_placeholder_png(path: Path, title: str) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_bytes(_MINIMAL_PNG_1X1)
  (path.with_suffix(".txt")).write_text(
    f"Placeholder figure: {title}\nmatplotlib unavailable: {_MPL_ERR}\n",
    encoding="utf-8",
  )


def _to_float(v: str) -> float | None:
  s = (v or "").strip()
  if not s:
    return None
  try:
    return float(s)
  except ValueError:
    return None


def _load_csv(path: Path) -> list[dict[str, str]]:
  if not path.exists():
    return []
  with path.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _plot_success_rate(summary_rows: list[dict[str, str]], out_path: Path) -> None:
  if plt is None:
    _write_placeholder_png(out_path, "Model Success Rate")
    return
  models = [r.get("model", "") for r in summary_rows]
  vals = [_to_float(r.get("success_rate", "")) or 0.0 for r in summary_rows]
  fig, ax = plt.subplots(figsize=(8, 4.5))
  if models:
    ax.bar(models, vals)
    ax.set_ylim(0.0, 1.0)
  else:
    ax.text(0.5, 0.5, "No data", ha="center", va="center")
  ax.set_title("Model Success Rate")
  ax.set_ylabel("success_rate")
  ax.set_xlabel("model")
  fig.tight_layout()
  fig.savefig(out_path, dpi=180)
  plt.close(fig)


def _plot_box(rows: list[dict[str, str]], field: str, title: str, y_label: str, out_path: Path) -> None:
  if plt is None:
    _write_placeholder_png(out_path, title)
    return
  by_model: dict[str, list[float]] = {}
  for r in rows:
    m = (r.get("model") or "").strip()
    x = _to_float(r.get(field, ""))
    if not m or x is None:
      continue
    by_model.setdefault(m, []).append(x)

  fig, ax = plt.subplots(figsize=(8, 4.5))
  models = sorted(by_model)
  data = [by_model[m] for m in models]
  if data:
    ax.boxplot(data, labels=models, showfliers=False)
  else:
    ax.text(0.5, 0.5, "No numeric data", ha="center", va="center")
  ax.set_title(title)
  ax.set_ylabel(y_label)
  ax.set_xlabel("model")
  fig.tight_layout()
  fig.savefig(out_path, dpi=180)
  plt.close(fig)


def _plot_heatmap(rows: list[dict[str, str]], out_path: Path) -> None:
  if plt is None:
    _write_placeholder_png(out_path, "Sample-Model Heatmap (CDR-H3 RMSD)")
    return
  models = sorted({(r.get("model") or "").strip() for r in rows if (r.get("model") or "").strip()})
  samples = sorted({(r.get("sample_id") or "").strip() for r in rows if (r.get("sample_id") or "").strip()})
  val_map: dict[tuple[str, str], float] = {}
  for r in rows:
    m = (r.get("model") or "").strip()
    s = (r.get("sample_id") or "").strip()
    x = _to_float(r.get("cdr_h3_rmsd", ""))
    if not m or not s or x is None:
      continue
    key = (s, m)
    if key not in val_map or x < val_map[key]:
      val_map[key] = x

  fig, ax = plt.subplots(figsize=(max(6, len(models) * 1.2), max(4, len(samples) * 0.2 + 2)))
  if models and samples and val_map:
    matrix: list[list[float]] = []
    for s in samples:
      row_vals = []
      for m in models:
        row_vals.append(val_map.get((s, m), math.nan))
      matrix.append(row_vals)
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(models, rotation=45, ha="right")
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples)
    ax.set_title("Sample-Model Heatmap (CDR-H3 RMSD)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
  else:
    ax.text(0.5, 0.5, "No RMSD data", ha="center", va="center")
    ax.set_title("Sample-Model Heatmap (CDR-H3 RMSD)")
    ax.set_xticks([])
    ax.set_yticks([])
  fig.tight_layout()
  fig.savefig(out_path, dpi=180)
  plt.close(fig)


def _plot_pareto(summary_rows: list[dict[str, str]], out_path: Path) -> None:
  if plt is None:
    _write_placeholder_png(out_path, "Pareto: Runtime vs CDR-H3 RMSD")
    return
  fig, ax = plt.subplots(figsize=(7.5, 5))
  has_points = False
  for r in summary_rows:
    model = (r.get("model") or "").strip()
    x = _to_float(r.get("median_runtime_sec", ""))
    y = _to_float(r.get("median_cdr_h3_rmsd", ""))
    success_rate = _to_float(r.get("success_rate", "")) or 0.0
    if not model or x is None or y is None:
      continue
    has_points = True
    size = 120 + 500 * success_rate
    ax.scatter([x], [y], s=size, alpha=0.7)
    ax.annotate(model, (x, y), xytext=(5, 5), textcoords="offset points")
  if not has_points:
    ax.text(0.5, 0.5, "No points (need runtime and RMSD)", ha="center", va="center")
  ax.set_title("Pareto: Runtime vs CDR-H3 RMSD")
  ax.set_xlabel("median_runtime_sec (lower is faster)")
  ax.set_ylabel("median_cdr_h3_rmsd (lower is better)")
  fig.tight_layout()
  fig.savefig(out_path, dpi=180)
  plt.close(fig)


def _write_report(summary_rows: list[dict[str, str]], report_path: Path) -> None:
  report_path.parent.mkdir(parents=True, exist_ok=True)
  lines: list[str] = []
  lines.append("# Evaluation Report (All Models)")
  lines.append("")
  lines.append("## Key Figures")
  lines.append("- `fig_01_success_rate.png`")
  lines.append("- `fig_02_cdr_h3_rmsd_boxplot.png`")
  lines.append("- `fig_03_runtime_boxplot.png`")
  lines.append("- `fig_04_sample_model_heatmap.png`")
  lines.append("- `fig_05_pareto_runtime_vs_rmsd.png`")
  lines.append("")
  lines.append("## Summary By Model")
  lines.append("")
  lines.append("| model | n_total | n_success | success_rate | median_runtime_sec | median_cdr_h3_rmsd |")
  lines.append("|---|---:|---:|---:|---:|---:|")
  for r in summary_rows:
    lines.append(
      "| {model} | {n_total} | {n_success} | {success_rate} | {median_runtime_sec} | {median_cdr_h3_rmsd} |".format(
        model=r.get("model", ""),
        n_total=r.get("n_total", ""),
        n_success=r.get("n_success", ""),
        success_rate=r.get("success_rate", ""),
        median_runtime_sec=r.get("median_runtime_sec", ""),
        median_cdr_h3_rmsd=r.get("median_cdr_h3_rmsd", ""),
      )
    )
  lines.append("")
  lines.append("## Notes")
  lines.append("- 失败样本已保留在评估长表中，不会被静默过滤。")
  lines.append("- 当结构结果缺失时，RMSD 相关图会显示空白/缺失。")
  if plt is None:
    lines.append(f"- 当前环境 matplotlib 不可用，图像已降级为占位图（原因：`{_MPL_ERR}`）。")
  report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
  parser = argparse.ArgumentParser(description="可视化 all-model 评估结果并生成报告")
  parser.add_argument(
    "--metrics-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/evaluation_long_metrics.csv"),
    help="包含 per-run 指标的长表",
  )
  parser.add_argument(
    "--summary-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/summary_by_model.csv"),
    help="模型汇总表",
  )
  parser.add_argument(
    "--fig-dir",
    type=Path,
    default=Path("outputs/evaluation/all_models/figures"),
    help="图像输出目录",
  )
  parser.add_argument(
    "--report-path",
    type=Path,
    default=Path("outputs/evaluation/all_models/report.md"),
    help="报告输出路径",
  )
  args = parser.parse_args()

  rows = _load_csv(args.metrics_csv)
  summary_rows = _load_csv(args.summary_csv)
  args.fig_dir.mkdir(parents=True, exist_ok=True)

  _plot_success_rate(summary_rows, args.fig_dir / "fig_01_success_rate.png")
  _plot_box(rows, "cdr_h3_rmsd", "CDR-H3 RMSD by Model", "cdr_h3_rmsd", args.fig_dir / "fig_02_cdr_h3_rmsd_boxplot.png")
  _plot_box(rows, "duration_sec", "Runtime by Model", "duration_sec", args.fig_dir / "fig_03_runtime_boxplot.png")
  _plot_heatmap(rows, args.fig_dir / "fig_04_sample_model_heatmap.png")
  _plot_pareto(summary_rows, args.fig_dir / "fig_05_pareto_runtime_vs_rmsd.png")
  _write_report(summary_rows, args.report_path)

  print(f"[OK] 图像输出目录: {args.fig_dir}")
  print(f"[OK] 报告: {args.report_path}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
