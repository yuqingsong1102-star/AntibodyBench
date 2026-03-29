#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


METRICS = [
  "cdr_h3_rmsd",
  "cdr_l3_rmsd",
  "rmsd",
  "tm_score",
  "dockq",
  "irmsd",
  "plddt_mean",
  "pae_mean",
  "runtime_sec",
  "gpu_mem_mb",
]


def _to_float(v: str) -> float | None:
  s = (v or "").strip()
  if not s:
    return None
  try:
    return float(s)
  except ValueError:
    return None


def aggregate(model: str, root: Path) -> Path:
  eval_dir = root / "outputs" / "evaluation" / model
  in_csv = eval_dir / "prediction_reference.csv"
  out_csv = eval_dir / "summary_metrics.csv"

  rows = []
  if in_csv.exists():
    with in_csv.open("r", encoding="utf-8", newline="") as f:
      rows = list(csv.DictReader(f))

  summary_rows: list[dict[str, str | int | float]] = []
  for metric in METRICS:
    values = []
    for r in rows:
      x = _to_float(r.get(metric, ""))
      if x is not None:
        values.append(x)
    if not values:
      summary_rows.append(
        {
          "metric": metric,
          "count": 0,
          "mean": "",
          "median": "",
          "stdev": "",
        }
      )
      continue
    summary_rows.append(
      {
        "metric": metric,
        "count": len(values),
        "mean": round(statistics.mean(values), 6),
        "median": round(statistics.median(values), 6),
        "stdev": round(statistics.pstdev(values), 6),
      }
    )

  eval_dir.mkdir(parents=True, exist_ok=True)
  with out_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["metric", "count", "mean", "median", "stdev"])
    writer.writeheader()
    writer.writerows(summary_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="聚合 prediction_reference.csv 的指标")
  parser.add_argument("--model", required=True, help="模型名称（如 germinal/boltzgen/BindCraft）")
  parser.add_argument(
    "--root",
    type=Path,
    default=Path(__file__).resolve().parents[1],
    help="AntibodyBench 根目录",
  )
  args = parser.parse_args()

  out_csv = aggregate(args.model, args.root)
  print(f"[OK] 写出指标汇总: {out_csv}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
