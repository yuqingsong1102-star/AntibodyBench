#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _has_atom_records(pdb_path: Path) -> bool:
  if not pdb_path.exists():
    return False
  with pdb_path.open("r", encoding="utf-8", errors="ignore") as f:
    for line in f:
      if line.startswith("ATOM") or line.startswith("HETATM"):
        return True
  return False


def build_index(interface_csv: Path, raw_pdb_dir: Path, out_csv: Path) -> tuple[int, int]:
  with interface_csv.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))

  out_rows = []
  skipped = 0
  for row in rows:
    pdb_id = (row.get("pdb_id") or row.get("pdb") or "").strip().lower()
    if not pdb_id:
      continue
    pdb_path = raw_pdb_dir / f"{pdb_id}.pdb"
    if not _has_atom_records(pdb_path):
      skipped += 1
      continue

    antibody_chain = (row.get("antibody_chain") or row.get("ab_chain") or row.get("ab_chains") or "").strip()
    antigen_chain = (row.get("antigen_chain") or row.get("ag_chain") or row.get("ag_chains") or "").strip()
    sample_id = f"{pdb_id}_{antibody_chain}_{antigen_chain}".strip("_")
    out_rows.append(
      {
        "sample_id": sample_id,
        "pdb_id": pdb_id,
        "reference_structure_path": str(pdb_path.resolve()),
        "antibody_chain": antibody_chain,
        "antigen_chain": antigen_chain,
      }
    )

  out_csv.parent.mkdir(parents=True, exist_ok=True)
  with out_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(
      f,
      fieldnames=["sample_id", "pdb_id", "reference_structure_path", "antibody_chain", "antigen_chain"],
    )
    writer.writeheader()
    writer.writerows(out_rows)
  return len(out_rows), skipped


def main() -> int:
  parser = argparse.ArgumentParser(description="构建抗体数据索引，并过滤空PDB")
  parser.add_argument("--interface-csv", required=True, type=Path)
  parser.add_argument("--raw-pdb-dir", required=True, type=Path)
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path(__file__).resolve().parents[1] / "inputs" / "antibody_datasets" / "dataset_index.csv",
  )
  args = parser.parse_args()

  kept, skipped = build_index(args.interface_csv, args.raw_pdb_dir, args.out_csv)
  print(f"[OK] 写出索引: {args.out_csv}")
  print(f"[INFO] 保留: {kept}, 过滤空PDB: {skipped}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
