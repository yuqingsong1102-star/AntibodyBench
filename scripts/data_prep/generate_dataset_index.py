#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
  parser = argparse.ArgumentParser(description="Generate AntibodyBench inputs/dataset_index.csv")
  parser.add_argument(
    "--raw-index",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "raw" / "dataset_index.csv",
    help="Raw dataset index (current repo format).",
  )
  parser.add_argument(
    "--reference-pdb-dir",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "raw" / "raw_pdbs",
    help="Directory containing reference PDBs named by base PDB id.",
  )
  parser.add_argument(
    "--out",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "inputs" / "antibody_datasets" / "dataset_index.csv",
    help="Output dataset index path compatible with AntibodyBench scripts/run.sh.",
  )
  args = parser.parse_args()

  args.out.parent.mkdir(parents=True, exist_ok=True)
  fieldnames = [
    "sample_id",
    "pdb_id",
    "reference_structure_path",
    "antigen_chain",
    "antibody_chain",
    "cdr_h3_start",
    "cdr_h3_end",
  ]

  with args.raw_index.open("r", encoding="utf-8", newline="") as f:
    raw_rows = list(csv.DictReader(f))

  out_rows: list[dict[str, str]] = []
  for r in raw_rows:
    raw_pdb_id = (r.get("pdb_id") or "").strip()
    antibody_chain = (r.get("antibody_chain") or "").strip()
    antigen_chain = (r.get("antigen_chain") or "").strip()
    if not raw_pdb_id or not antibody_chain or not antigen_chain:
      continue

    pdb_base = raw_pdb_id.split("-")[0].strip().lower()
    ref_path = args.reference_pdb_dir / f"{pdb_base}.pdb"
    sample_id = f"{pdb_base}_{antibody_chain}_{antigen_chain}"
    out_rows.append(
      {
        "sample_id": sample_id,
        "pdb_id": pdb_base,
        "reference_structure_path": str(ref_path),
        "antigen_chain": antigen_chain,
        "antibody_chain": antibody_chain,
        "cdr_h3_start": "",
        "cdr_h3_end": "",
      }
    )

  with args.out.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(out_rows)
  print(f"[OK] Wrote dataset index: {args.out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

