#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import csv
import io
from pathlib import Path

_MPL_ERR = ""
try:
  with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore
except Exception as exc:  # pragma: no cover
  plt = None
  _MPL_ERR = f"{type(exc).__name__}: {exc}"


_MINIMAL_PNG_1X1 = (
  b"\x89PNG\r\n\x1a\n"
  b"\x00\x00\x00\rIHDR"
  b"\x00\x00\x00\x01\x00\x00\x00\x01\x08\x06\x00\x00\x00"
  b"\x1f\x15\xc4\x89"
  b"\x00\x00\x00\nIDATx\x9cc`\x00\x00\x00\x02\x00\x01"
  b"\xe2!\xbc3"
  b"\x00\x00\x00\x00IEND\xaeB`\x82"
)


def _load_csv(path: Path) -> list[dict[str, str]]:
  if not path.exists():
    return []
  with path.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _to_float(value: str) -> float | None:
  text = (value or "").strip()
  if not text:
    return None
  try:
    return float(text)
  except ValueError:
    return None


def _write_placeholder_png(path: Path, title: str) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_bytes(_MINIMAL_PNG_1X1)
  (path.with_suffix(".txt")).write_text(
    f"Placeholder figure: {title}\nmatplotlib unavailable: {_MPL_ERR}\n",
    encoding="utf-8",
  )


def _plot_bar(summary_rows: list[dict[str, str]], field: str, title: str, ylabel: str, out_path: Path) -> None:
  if plt is None:
    _write_placeholder_png(out_path, title)
    return
  models = [row.get("model", "") for row in summary_rows]
  values = [_to_float(row.get(field, "")) or 0.0 for row in summary_rows]
  fig, ax = plt.subplots(figsize=(8, 4.5))
  if models:
    ax.bar(models, values)
  else:
    ax.text(0.5, 0.5, "No data", ha="center", va="center")
  ax.set_title(title)
  ax.set_ylabel(ylabel)
  ax.set_xlabel("model")
  fig.tight_layout()
  fig.savefig(out_path, dpi=180)
  plt.close(fig)


def _plot_pass_at_k(summary_rows: list[dict[str, str]], out_path: Path) -> None:
  if plt is None:
    _write_placeholder_png(out_path, "Bindability Pass@K")
    return
  models = [row.get("model", "") for row in summary_rows]
  x = list(range(len(models)))
  width = 0.22
  vals1 = [_to_float(row.get("bindability_pass_at_1", "")) or 0.0 for row in summary_rows]
  vals5 = [_to_float(row.get("bindability_pass_at_5", "")) or 0.0 for row in summary_rows]
  vals10 = [_to_float(row.get("bindability_pass_at_10", "")) or 0.0 for row in summary_rows]
  fig, ax = plt.subplots(figsize=(9, 4.8))
  ax.bar([idx - width for idx in x], vals1, width=width, label="Pass@1")
  ax.bar(x, vals5, width=width, label="Pass@5")
  ax.bar([idx + width for idx in x], vals10, width=width, label="Pass@10")
  ax.set_xticks(x)
  ax.set_xticklabels(models)
  ax.set_ylim(0.0, 1.0)
  ax.set_ylabel("rate")
  ax.set_title("Bindability Pass@K")
  ax.legend()
  fig.tight_layout()
  fig.savefig(out_path, dpi=180)
  plt.close(fig)


def _plot_scatter(sequence_rows: list[dict[str, str]], x_field: str, y_field: str, title: str, xlabel: str, ylabel: str, out_path: Path) -> None:
  if plt is None:
    _write_placeholder_png(out_path, title)
    return
  xs = []
  ys = []
  for row in sequence_rows:
    x = _to_float(row.get(x_field, ""))
    y = _to_float(row.get(y_field, ""))
    if x is None or y is None:
      continue
    xs.append(x)
    ys.append(y)
  fig, ax = plt.subplots(figsize=(6.8, 5.2))
  if xs:
    ax.scatter(xs, ys, alpha=0.7)
  else:
    ax.text(0.5, 0.5, "No numeric data", ha="center", va="center")
  ax.set_title(title)
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  fig.tight_layout()
  fig.savefig(out_path, dpi=180)
  plt.close(fig)


def _plot_box(sequence_rows: list[dict[str, str]], field: str, title: str, ylabel: str, out_path: Path) -> None:
  if plt is None:
    _write_placeholder_png(out_path, title)
    return
  by_model: dict[str, list[float]] = {}
  for row in sequence_rows:
    model = (row.get("model") or "").strip()
    value = _to_float(row.get(field, ""))
    if not model or value is None:
      continue
    by_model.setdefault(model, []).append(value)
  fig, ax = plt.subplots(figsize=(8, 4.5))
  models = sorted(by_model)
  data = [by_model[model] for model in models]
  if data:
    ax.boxplot(data, labels=models, showfliers=False)
  else:
    ax.text(0.5, 0.5, "No numeric data", ha="center", va="center")
  ax.set_title(title)
  ax.set_ylabel(ylabel)
  ax.set_xlabel("model")
  fig.tight_layout()
  fig.savefig(out_path, dpi=180)
  plt.close(fig)


def generate_external_figures(sequence_csv: Path, model_summary_csv: Path, fig_dir: Path) -> Path:
  sequence_rows = _load_csv(sequence_csv)
  summary_rows = _load_csv(model_summary_csv)
  fig_dir.mkdir(parents=True, exist_ok=True)
  _plot_bar(summary_rows, "viability_rate", "Viability Rate", "viability_rate", fig_dir / "fig_01_viability_rate.png")
  _plot_pass_at_k(summary_rows, fig_dir / "fig_02_bindability_pass_at_k.png")
  _plot_box(sequence_rows, "median_pdockq2", "pDockQ2 by Model", "pdockq2", fig_dir / "fig_03_pdockq2_boxplot.png")
  _plot_box(sequence_rows, "median_external_iptm", "External iPTM by Model", "external_iptm", fig_dir / "fig_04_external_iptm_boxplot.png")
  _plot_scatter(
    sequence_rows,
    "median_interface_bsa",
    "median_interface_dG",
    "Interface BSA vs dG",
    "interface_bsa",
    "interface_dG",
    fig_dir / "fig_05_interface_bsa_dg_scatter.png",
  )
  _plot_scatter(
    sequence_rows,
    "median_surface_hydrophobicity",
    "median_binder_score",
    "Surface Hydrophobicity vs Binder Score",
    "surface_hydrophobicity",
    "binder_score",
    fig_dir / "fig_06_surface_hydrophobicity_binder_score_scatter.png",
  )
  _plot_bar(summary_rows, "robust_pass_rate", "Robust Pass Rate", "robust_pass_rate", fig_dir / "fig_07_robust_pass_rate.png")
  return fig_dir


def main() -> int:
  parser = argparse.ArgumentParser(description="绘制 Phase1 外部评估图表")
  parser.add_argument(
    "--sequence-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_sequence_level.csv"),
    help="candidate 汇总表",
  )
  parser.add_argument(
    "--model-summary-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_model_summary.csv"),
    help="模型汇总表",
  )
  parser.add_argument(
    "--fig-dir",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/figures"),
    help="图表输出目录",
  )
  args = parser.parse_args()
  out_dir = generate_external_figures(args.sequence_csv, args.model_summary_csv, args.fig_dir)
  print(f"[OK] 图表目录: {out_dir}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())