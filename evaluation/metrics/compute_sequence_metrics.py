#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from evaluation.schema import EVALUATION_LONG_FIELDS


AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def _to_int(v: str) -> int | None:
  s = (v or "").strip()
  if not s:
    return None
  try:
    return int(float(s))
  except ValueError:
    return None


def _percentile(values: list[float], q: float) -> float | None:
  if not values:
    return None
  vals = sorted(values)
  if len(vals) == 1:
    return vals[0]
  pos = (len(vals) - 1) * q
  low = int(pos)
  high = min(low + 1, len(vals) - 1)
  if low == high:
    return vals[low]
  frac = pos - low
  return vals[low] * (1.0 - frac) + vals[high] * frac


def _round_or_blank(x: float | None) -> str:
  if x is None:
    return ""
  return str(round(x, 6))


def _compute_seq_fields(row: dict[str, str]) -> dict[str, str]:
  seq = (row.get("pred_sequence") or "").strip().upper()
  if not seq:
    row["pred_seq_len"] = ""
    row["aa_valid_ratio"] = ""
    row["cys_count"] = ""
    row["motif_h3_len_in_8_20"] = ""
    row["has_sequence"] = "0"
    return row

  n = len(seq)
  valid_n = sum(1 for c in seq if c in AA20)
  cys_count = seq.count("C")
  row["pred_seq_len"] = str(n)
  row["aa_valid_ratio"] = _round_or_blank(valid_n / n if n > 0 else None)
  row["cys_count"] = str(cys_count)
  row["motif_h3_len_in_8_20"] = "1" if 8 <= n <= 20 else "0"
  row["has_sequence"] = "1"
  return row


def _load_rows(in_csv: Path) -> list[dict[str, str]]:
  if not in_csv.exists():
    return []
  with in_csv.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)


def compute_sequence_metrics(
  in_csv: Path,
  out_csv: Path,
) -> Path:
  rows = _load_rows(in_csv)
  out_rows: list[dict[str, str]] = []
  for row in rows:
    base = {k: (row.get(k, "") or "") for k in EVALUATION_LONG_FIELDS}
    base = _compute_seq_fields(base)
    out_rows.append(base)

  _write_csv(out_csv, EVALUATION_LONG_FIELDS, out_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="计算 Level-0/1 序列指标并回写长表")
  parser.add_argument(
    "--in-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/evaluation_long.csv"),
    help="输入长表路径",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/evaluation_long_seq.csv"),
    help="输出带序列指标的长表路径",
  )
  args = parser.parse_args()

  out_csv = compute_sequence_metrics(in_csv=args.in_csv, out_csv=args.out_csv)
  print(f"[OK] 序列指标长表: {out_csv}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

