#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path

from schema import SUMMARY_BY_MODEL_FIELDS


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


def _round_or_blank(v: float | None) -> str:
  if v is None:
    return ""
  return str(round(v, 6))


def _legacy_aggregate(model: str, root: Path) -> Path:
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


def _aggregate_from_long(in_csv: Path, out_csv: Path) -> Path:
  rows: list[dict[str, str]] = []
  if in_csv.exists():
    with in_csv.open("r", encoding="utf-8", newline="") as f:
      rows = list(csv.DictReader(f))

  by_model: dict[str, list[dict[str, str]]] = {}
  for r in rows:
    model = (r.get("model") or "").strip() or "unknown"
    by_model.setdefault(model, []).append(r)

  summary_rows: list[dict[str, str]] = []
  for model in sorted(by_model):
    rs = by_model[model]
    n_total = len(rs)
    n_success = sum(1 for r in rs if (r.get("status") or "").strip() == "ok")
    has_structure = sum(1 for r in rs if (r.get("has_structure") or "").strip() == "1")
    has_sequence = sum(1 for r in rs if (r.get("has_sequence") or "").strip() == "1")
    runtime_values = [_to_float(r.get("duration_sec", "")) for r in rs]
    runtime_values = [x for x in runtime_values if x is not None]
    seq_len_values = [_to_float(r.get("pred_seq_len", "")) for r in rs]
    seq_len_values = [x for x in seq_len_values if x is not None]
    rmsd_values = [_to_float(r.get("cdr_h3_rmsd", "")) for r in rs]
    rmsd_values = [x for x in rmsd_values if x is not None]
    rmsd_values.sort()

    row = {k: "" for k in SUMMARY_BY_MODEL_FIELDS}
    row["model"] = model
    row["n_total"] = str(n_total)
    row["n_success"] = str(n_success)
    row["success_rate"] = _round_or_blank((n_success / n_total) if n_total else None)
    row["has_structure_rate"] = _round_or_blank((has_structure / n_total) if n_total else None)
    row["has_sequence_rate"] = _round_or_blank((has_sequence / n_total) if n_total else None)
    row["median_runtime_sec"] = _round_or_blank(statistics.median(runtime_values) if runtime_values else None)
    row["median_seq_len"] = _round_or_blank(statistics.median(seq_len_values) if seq_len_values else None)
    row["median_cdr_h3_rmsd"] = _round_or_blank(statistics.median(rmsd_values) if rmsd_values else None)
    if rmsd_values:
      p25_idx = int(0.25 * (len(rmsd_values) - 1))
      p75_idx = int(0.75 * (len(rmsd_values) - 1))
      row["p25_cdr_h3_rmsd"] = _round_or_blank(rmsd_values[p25_idx])
      row["p75_cdr_h3_rmsd"] = _round_or_blank(rmsd_values[p75_idx])
    summary_rows.append(row)

  out_csv.parent.mkdir(parents=True, exist_ok=True)
  with out_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=SUMMARY_BY_MODEL_FIELDS)
    writer.writeheader()
    writer.writerows(summary_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="聚合评估指标（支持 long table 与 legacy 两种模式）")
  parser.add_argument("--model", help="legacy 模式下模型名称（如 germinal/boltzgen/BindCraft）")
  parser.add_argument(
    "--root",
    type=Path,
    default=Path(__file__).resolve().parents[2],
    help="AntibodyBench 根目录",
  )
  parser.add_argument(
    "--in-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/evaluation_long_metrics.csv"),
    help="long table 输入路径",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/summary_by_model.csv"),
    help="long table 模式输出路径",
  )
  parser.add_argument(
    "--mode",
    choices=["long", "legacy"],
    default="long",
    help="聚合模式",
  )
  args = parser.parse_args()

  if args.mode == "legacy":
    if not args.model:
      raise SystemExit("--mode legacy 时必须提供 --model")
    out_csv = _legacy_aggregate(args.model, args.root)
    print(f"[OK] 写出 legacy 指标汇总: {out_csv}")
    return 0

  out_csv = _aggregate_from_long(args.in_csv, args.out_csv)
  print(f"[OK] 写出 long table 模型汇总: {out_csv}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
