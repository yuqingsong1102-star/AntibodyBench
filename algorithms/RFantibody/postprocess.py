#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
  parser = argparse.ArgumentParser(description="RFantibody 后处理")
  parser.add_argument("--dataset-csv", required=True, type=Path)
  parser.add_argument("--prediction-dir", required=True, type=Path)
  parser.add_argument("--out-dir", required=True, type=Path)
  args = parser.parse_args()

  args.out_dir.mkdir(parents=True, exist_ok=True)
  rows = []
  with args.dataset_csv.open("r", encoding="utf-8", newline="") as f:
    for row in csv.DictReader(f):
      sample_id = (row.get("sample_id") or "").strip()
      if not sample_id:
        continue
      pred = args.prediction_dir / f"{sample_id}.pdb"
      if not pred.exists():
        continue
      rows.append(
        {
          "pdb_id": row.get("pdb_id", ""),
          "seed": 0,
          "sample": sample_id,
          "ranking_score": "",
          "prediction_path": str(pred.resolve()),
          "cdr_h3_rmsd": "",
          "cdr_l3_rmsd": "",
          "rmsd": "",
          "tm_score": "",
          "dockq": "",
          "irmsd": "",
          "plddt_mean": "",
          "pae_mean": "",
          "runtime_sec": "",
          "gpu_mem_mb": "",
        }
      )

  out_csv = args.out_dir / "prediction_reference.csv"
  with out_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(
      f,
      fieldnames=[
        "pdb_id",
        "seed",
        "sample",
        "ranking_score",
        "prediction_path",
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
      ],
    )
    writer.writeheader()
    writer.writerows(rows)
  print(f"[OK] RFantibody postprocess 输出: {out_csv}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
