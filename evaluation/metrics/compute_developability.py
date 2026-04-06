#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

from evaluation.schema import EXTERNAL_SEED_FIELDS


GLYCO_RE = re.compile(r"N[^P][ST]")
DEAMIDATION_RE = re.compile(r"N[GS]")
ISOMERIZATION_RE = re.compile(r"D[GST]")

_AA_CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1}


def _load_csv(path: Path) -> list[dict[str, str]]:
  if not path.exists():
    return []
  with path.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=EXTERNAL_SEED_FIELDS)
    writer.writeheader()
    writer.writerows(rows)


def _round_or_blank(value: float | None) -> str:
  if value is None:
    return ""
  return str(round(value, 6))


def _motif_count(pattern: re.Pattern[str], seq: str) -> int:
  return sum(1 for _ in pattern.finditer(seq))


def _sequence_entropy(seq: str) -> float | None:
  if not seq:
    return None
  counts = {}
  for aa in seq:
    counts[aa] = counts.get(aa, 0) + 1
  total = len(seq)
  entropy = 0.0
  for count in counts.values():
    p = count / total
    entropy -= p * math.log(p, 2)
  return entropy


def _net_charge(seq: str) -> float | None:
  if not seq:
    return None
  return float(sum(_AA_CHARGE.get(aa, 0.0) for aa in seq))


def _score_row(row: dict[str, str]) -> dict[str, str]:
  out = {k: (row.get(k, "") or "") for k in EXTERNAL_SEED_FIELDS}
  seq = (out.get("pred_sequence") or "").strip().upper()

  cys_count = int((out.get("cys_count") or "0") or 0) if seq else 0
  free_cys_flag = cys_count > 0
  glyco_count = _motif_count(GLYCO_RE, seq)
  deamidation_count = _motif_count(DEAMIDATION_RE, seq)
  isomerization_count = _motif_count(ISOMERIZATION_RE, seq)
  liability_motif_count = glyco_count + deamidation_count + isomerization_count

  out["free_cys_flag"] = "1" if free_cys_flag else "0"
  out["glycosylation_motif_flag"] = "1" if glyco_count > 0 else "0"
  out["deamidation_risk_flag"] = "1" if deamidation_count > 0 else "0"
  out["isomerization_risk_flag"] = "1" if isomerization_count > 0 else "0"
  out["liability_motif_count"] = str(liability_motif_count)
  out["sequence_entropy"] = _round_or_blank(_sequence_entropy(seq))
  out["net_charge"] = _round_or_blank(_net_charge(seq))
  basic_liability = free_cys_flag or liability_motif_count > 0
  out["basic_liability_flag"] = "1" if basic_liability else "0"

  surface_hydrophobicity = None
  if out.get("surface_hydrophobicity"):
    try:
      surface_hydrophobicity = float(out["surface_hydrophobicity"])
    except Exception:
      surface_hydrophobicity = None

  checks = [not basic_liability]
  if surface_hydrophobicity is not None:
    checks.append(surface_hydrophobicity <= 0.45)
  out["developability_pass"] = "1" if all(checks) else "0"
  return out


def compute_developability(in_csv: Path, out_csv: Path) -> Path:
  rows = _load_csv(in_csv)
  out_rows = [_score_row(row) for row in rows]
  _write_csv(out_csv, out_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="计算 Developability 与 Novelty 相关字段")
  parser.add_argument(
    "--in-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="seed 级主表",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="输出 seed 级主表（可覆盖输入）",
  )
  args = parser.parse_args()
  out = compute_developability(args.in_csv, args.out_csv)
  print(f"[OK] developability: {out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
