#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> int:
  parser = argparse.ArgumentParser(description="template_model 预处理")
  parser.add_argument("--af3-input", required=True, type=Path)
  parser.add_argument("--dataset-csv", required=True, type=Path)
  parser.add_argument("--out-dir", required=True, type=Path)
  args = parser.parse_args()

  args.out_dir.mkdir(parents=True, exist_ok=True)
  with args.dataset_csv.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))
  for row in rows:
    sample_id = (row.get("sample_id") or "").strip()
    if not sample_id:
      continue
    payload = {
      "sample_id": sample_id,
      "pdb_id": (row.get("pdb_id") or "").strip(),
      "reference_structure_path": (row.get("reference_structure_path") or "").strip(),
    }
    (args.out_dir / f"{sample_id}.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
  print(f"[OK] template_model preprocess 输出: {args.out_dir}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
