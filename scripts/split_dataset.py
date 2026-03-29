#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path


def split_dataset(index_csv: Path, out_dir: Path, test_ratio: float, seed: int) -> tuple[Path, Path]:
  with index_csv.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))
  if not rows:
    out_dir.mkdir(parents=True, exist_ok=True)
    train_csv = out_dir / "train.csv"
    test_csv = out_dir / "test.csv"
    for p in (train_csv, test_csv):
      with p.open("w", encoding="utf-8", newline="") as wf:
        writer = csv.writer(wf)
        writer.writerow(["sample_id", "pdb_id", "reference_structure_path", "antibody_chain", "antigen_chain"])
    return train_csv, test_csv

  by_pdb: dict[str, list[dict[str, str]]] = {}
  for r in rows:
    pdb_id = (r.get("pdb_id") or "").strip().lower()
    if not pdb_id:
      continue
    by_pdb.setdefault(pdb_id, []).append(r)

  pdb_ids = sorted(by_pdb.keys())
  random.Random(seed).shuffle(pdb_ids)
  test_n = max(1, int(len(pdb_ids) * test_ratio)) if pdb_ids else 0
  test_set = set(pdb_ids[:test_n])

  train_rows: list[dict[str, str]] = []
  test_rows: list[dict[str, str]] = []
  for pdb_id, group in by_pdb.items():
    (test_rows if pdb_id in test_set else train_rows).extend(group)

  out_dir.mkdir(parents=True, exist_ok=True)
  train_csv = out_dir / "train.csv"
  test_csv = out_dir / "test.csv"
  fieldnames = list(rows[0].keys())
  for p, data in ((train_csv, train_rows), (test_csv, test_rows)):
    with p.open("w", encoding="utf-8", newline="") as wf:
      writer = csv.DictWriter(wf, fieldnames=fieldnames)
      writer.writeheader()
      writer.writerows(data)
  return train_csv, test_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="按 pdb_id 划分 train/test，避免泄漏")
  parser.add_argument("--index-csv", required=True, type=Path)
  parser.add_argument("--out-dir", required=True, type=Path)
  parser.add_argument("--test-ratio", type=float, default=0.2)
  parser.add_argument("--seed", type=int, default=42)
  args = parser.parse_args()
  train_csv, test_csv = split_dataset(args.index_csv, args.out_dir, args.test_ratio, args.seed)
  print(f"[OK] train: {train_csv}")
  print(f"[OK] test : {test_csv}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
