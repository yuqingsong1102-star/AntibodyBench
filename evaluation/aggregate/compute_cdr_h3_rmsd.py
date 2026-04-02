#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
  sys.path.insert(0, str(ROOT))

from evaluation.cdr_rmsd import compute_antigen_aligned_cdr_rmsd  # noqa: E402
from schema import EVALUATION_LONG_FIELDS  # noqa: E402


def _to_int(v: str) -> int | None:
  s = (v or "").strip()
  if not s:
    return None
  try:
    return int(float(s))
  except ValueError:
    return None


def _load_rows(in_csv: Path) -> list[dict[str, str]]:
  if not in_csv.exists():
    return []
  with in_csv.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _write_rows(out_csv: Path, rows: list[dict[str, str]]) -> None:
  out_csv.parent.mkdir(parents=True, exist_ok=True)
  with out_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=EVALUATION_LONG_FIELDS)
    writer.writeheader()
    writer.writerows(rows)


def _is_eligible(row: dict[str, str]) -> tuple[bool, str]:
  structure_path = Path((row.get("structure_path") or "").strip()) if (row.get("structure_path") or "").strip() else None
  ref_path = Path((row.get("reference_complex_path") or "").strip()) if (row.get("reference_complex_path") or "").strip() else None
  antigen_chain = (row.get("antigen_chain_id") or "").strip()
  antibody_chain = (row.get("antibody_chain_id") or "").strip()
  cdr_status = (row.get("cdr_h3_status") or "").strip().lower()
  cdr_start = _to_int(row.get("cdr_h3_start", ""))
  cdr_end = _to_int(row.get("cdr_h3_end", ""))

  if not structure_path or not structure_path.exists():
    return False, "missing_pred_structure"
  if not ref_path or not ref_path.exists():
    return False, "missing_reference_complex"
  if cdr_status not in {"ok", ""}:
    return False, f"cdr_h3_status={cdr_status}"
  if cdr_start is None or cdr_end is None:
    return False, "missing_cdr_h3_range"
  if cdr_start > cdr_end:
    return False, "invalid_cdr_h3_range"
  if not antigen_chain or not antibody_chain:
    return False, "missing_chain_id"
  return True, ""


def compute_cdr_h3_rmsd(in_csv: Path, out_csv: Path) -> Path:
  rows = _load_rows(in_csv)
  out_rows: list[dict[str, str]] = []

  for row in rows:
    cur = {k: (row.get(k, "") or "") for k in EVALUATION_LONG_FIELDS}
    ok, reason = _is_eligible(cur)
    cur["struct_eval_eligible"] = "1" if ok else "0"
    cur["metric_error"] = ""
    cur["cdr_h3_rmsd"] = ""
    cur["cdr_h3_atom_count"] = ""
    if not ok:
      cur["metric_error"] = reason
      out_rows.append(cur)
      continue

    try:
      result = compute_antigen_aligned_cdr_rmsd(
        pred_path=Path(cur["structure_path"]),
        ref_path=Path(cur["reference_complex_path"]),
        antigen_chain_id=cur["antigen_chain_id"],
        antibody_chain_id=cur["antibody_chain_id"],
        cdr_start=int(float(cur["cdr_h3_start"])),
        cdr_end=int(float(cur["cdr_h3_end"])),
      )
      if result is None:
        cur["metric_error"] = "rmsd_result_none"
      else:
        cur["cdr_h3_rmsd"] = str(round(float(result.rmsd), 6))
        cur["cdr_h3_atom_count"] = str(int(result.atom_count))
    except Exception as e:
      cur["metric_error"] = f"rmsd_exception:{type(e).__name__}:{e}"
    out_rows.append(cur)

  _write_rows(out_csv, out_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="批量计算 cdr_h3_rmsd 并回写长表")
  parser.add_argument(
    "--in-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/evaluation_long_seq.csv"),
    help="输入长表路径（建议先跑序列指标脚本）",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/evaluation_long_metrics.csv"),
    help="输出路径",
  )
  args = parser.parse_args()

  out_csv = compute_cdr_h3_rmsd(in_csv=args.in_csv, out_csv=args.out_csv)
  print(f"[OK] 写出结构指标长表: {out_csv}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

