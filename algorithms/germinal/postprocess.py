#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


class PostProcess:
  def __init__(self, dataset_csv: Path, prediction_dir: Path, out_dir: Path) -> None:
    self.dataset_csv = dataset_csv
    self.prediction_dir = prediction_dir
    self.out_dir = out_dir

  def _extract_cdr_region(self, src_struct: Path, out_struct: Path) -> None:
    out_struct.write_text(src_struct.read_text(encoding="utf-8"), encoding="utf-8")

  def postprocess(self) -> None:
    self.out_dir.mkdir(parents=True, exist_ok=True)
    cdr_dir = self.out_dir / "cdr_regions"
    cdr_dir.mkdir(parents=True, exist_ok=True)

    with self.dataset_csv.open("r", encoding="utf-8", newline="") as f:
      dataset_rows = list(csv.DictReader(f))

    rows = []
    for row in dataset_rows:
      sample_id = str(row.get("sample_id", "")).strip()
      if not sample_id:
        continue
      pred_path = None
      for ext in ("pdb", "cif"):
        cand = self.prediction_dir / f"{sample_id}.{ext}"
        if cand.exists():
          pred_path = cand
          break
      if pred_path is None:
        continue

      cdr_h3 = cdr_dir / f"{sample_id}_cdr_h3.{pred_path.suffix.lstrip('.')}"
      cdr_l3 = cdr_dir / f"{sample_id}_cdr_l3.{pred_path.suffix.lstrip('.')}"
      self._extract_cdr_region(pred_path, cdr_h3)
      self._extract_cdr_region(pred_path, cdr_l3)

      rows.append(
        {
          "pdb_id": row.get("pdb_id", ""),
          "seed": 0,
          "sample": sample_id,
          "ranking_score": "",
          "prediction_path": str(pred_path.resolve()),
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

    out_csv = self.out_dir / "prediction_reference.csv"
    with out_csv.open("w", encoding="utf-8", newline="") as f:
      fieldnames = [
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
      ]
      writer = csv.DictWriter(f, fieldnames=fieldnames)
      writer.writeheader()
      writer.writerows(rows)
    print(f"[OK] germinal postprocess 输出: {out_csv}")


def main() -> int:
  parser = argparse.ArgumentParser(description="germinal 后处理")
  parser.add_argument("--dataset-csv", required=True, type=Path)
  parser.add_argument("--prediction-dir", required=True, type=Path)
  parser.add_argument("--out-dir", required=True, type=Path)
  args = parser.parse_args()
  PostProcess(args.dataset_csv, args.prediction_dir, args.out_dir).postprocess()
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
