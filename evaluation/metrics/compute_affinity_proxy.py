#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from evaluation.schema import EXTERNAL_SEED_FIELDS


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


def _to_float(text: str) -> float | None:
  text = (text or "").strip()
  if not text:
    return None
  try:
    return float(text)
  except Exception:
    return None


def _round_or_blank(value: float | None) -> str:
  if value is None:
    return ""
  return str(round(value, 6))


def _score_row(row: dict[str, str]) -> dict[str, str]:
  out = {k: (row.get(k, "") or "") for k in EXTERNAL_SEED_FIELDS}
  if out.get("viability_pass") != "1" or out.get("binding_plausibility_pass") != "1":
    out["affinity_proxy_pass"] = "0"
    return out

  contact_pairs = _to_float(out.get("interface_contact_pair_count", ""))
  interface_bsa = _to_float(out.get("interface_bsa", ""))
  clash_ratio = _to_float(out.get("interface_clash_ratio", ""))
  binder_fraction = _to_float(out.get("interface_binder_fraction", ""))
  pdockq2 = _to_float(out.get("pdockq2", ""))
  surf_hydro = _to_float(out.get("surface_hydrophobicity", ""))

  if contact_pairs is not None:
    hbond_count = max(0, int(round(contact_pairs * 0.08)))
    saltbridge_count = max(0, int(round(contact_pairs * 0.02)))
  else:
    hbond_count = 0
    saltbridge_count = 0
  unsat_hbond = max(0, int(round(hbond_count * 0.2)))

  interface_dSASA = interface_bsa
  interface_dG = None
  if interface_bsa is not None:
    clash_penalty = (clash_ratio or 0.0) * 50.0
    interface_dG = (-0.018 * interface_bsa) - (0.25 * hbond_count) - (0.4 * saltbridge_count) + clash_penalty
  interface_dG_dSASA_ratio = (interface_dG / interface_dSASA) if (interface_dG is not None and interface_dSASA and interface_dSASA > 0) else None

  interface_hydrophobicity = None
  if binder_fraction is not None:
    interface_hydrophobicity = min(1.0, max(0.0, binder_fraction / 100.0))

  binder_score = None
  if pdockq2 is not None and interface_dG is not None and interface_bsa is not None:
    binder_score = (pdockq2 * 100.0) + max(0.0, -interface_dG) + (interface_bsa / 50.0)

  sap_score = None
  if surf_hydro is not None:
    sap_score = surf_hydro * 100.0

  out["interface_hbond_count"] = str(hbond_count)
  out["interface_saltbridge_count"] = str(saltbridge_count)
  out["interface_unsat_hbond_count"] = str(unsat_hbond)
  out["interface_dSASA"] = _round_or_blank(interface_dSASA)
  out["interface_dG"] = _round_or_blank(interface_dG)
  out["interface_dG_dSASA_ratio"] = _round_or_blank(interface_dG_dSASA_ratio)
  out["interface_hydrophobicity"] = _round_or_blank(interface_hydrophobicity)
  out["binder_score"] = _round_or_blank(binder_score)
  out["sap_score"] = _round_or_blank(sap_score)

  affinity_checks = []
  if interface_dG is not None:
    affinity_checks.append(interface_dG <= -5.0)
  if interface_dG_dSASA_ratio is not None:
    affinity_checks.append(interface_dG_dSASA_ratio <= -0.01)
  affinity_checks.append(hbond_count >= 2)
  if surf_hydro is not None:
    affinity_checks.append(surf_hydro <= 0.45)
  out["affinity_proxy_pass"] = "1" if (affinity_checks and all(affinity_checks)) else "0"
  return out


def compute_affinity_proxy(in_csv: Path, out_csv: Path) -> Path:
  rows = _load_csv(in_csv)
  out_rows = [_score_row(row) for row in rows]
  _write_csv(out_csv, out_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="计算 Affinity Proxy 评分层")
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
  out = compute_affinity_proxy(args.in_csv, args.out_csv)
  print(f"[OK] affinity proxy: {out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
