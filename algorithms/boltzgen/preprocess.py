#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


class PreProcess:
  def __init__(self, af3_input: Path, dataset_csv: Path, out_dir: Path) -> None:
    self.af3_input = af3_input
    self.dataset_csv = dataset_csv
    self.out_dir = out_dir

  def _load_af3(self) -> dict[str, dict[str, Any]]:
    if not self.af3_input.exists():
      return {}
    with self.af3_input.open("r", encoding="utf-8") as f:
      data = json.load(f)
    if not isinstance(data, list):
      return {}
    out: dict[str, dict[str, Any]] = {}
    for item in data:
      if isinstance(item, dict) and item.get("sample_id"):
        out[str(item["sample_id"]).strip()] = item
    return out

  def preprocess(self) -> None:
    self.out_dir.mkdir(parents=True, exist_ok=True)
    af3_map = self._load_af3()
    with self.dataset_csv.open("r", encoding="utf-8", newline="") as f:
      rows = list(csv.DictReader(f))

    manifest_records: list[dict[str, Any]] = []
    skip_reasons: dict[str, int] = {}
    ready_count = 0

    for row in rows:
      sample_id = str(row.get("sample_id", "")).strip()
      pdb_id = str(row.get("pdb_id", "")).strip().lower()
      ref_path_raw = str(row.get("reference_structure_path", "")).strip()
      ref_path = Path(ref_path_raw) if ref_path_raw else None

      skip_reason = ""
      if not sample_id:
        skip_reason = "missing_sample_id"
      elif not pdb_id:
        skip_reason = "missing_pdb_id"
      elif not ref_path_raw:
        skip_reason = "missing_reference_structure_path"
      elif ref_path is not None and not ref_path.exists():
        skip_reason = "reference_structure_not_found"

      if skip_reason:
        skip_reasons[skip_reason] = skip_reasons.get(skip_reason, 0) + 1
        manifest_records.append(
          {
            "sample_id": sample_id,
            "pdb_id": pdb_id,
            "input_json": "",
            "status": "skip",
            "reason": skip_reason,
          }
        )
        continue

      af3 = af3_map.get(sample_id, {})
      input_json_path = self.out_dir / f"{sample_id}.json"
      payload = {
        "sample_id": sample_id,
        "pdb_id": pdb_id,
        "reference_structure_path": ref_path_raw,
        "antibody_chain": str(row.get("antibody_chain", "")).strip(),
        "antigen_chain": str(row.get("antigen_chain", "")).strip(),
        "cdr": af3.get("cdr", {}),
        "heavy_sequence": af3.get("heavy_sequence", ""),
        "light_sequence": af3.get("light_sequence", ""),
      }
      input_json_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
      )
      (self.out_dir / f"{sample_id}.a3m").write_text(
        ">query\n" + (payload["heavy_sequence"] or "") + "\n",
        encoding="utf-8",
      )
      ready_count += 1
      manifest_records.append(
        {
          "sample_id": sample_id,
          "pdb_id": pdb_id,
          "input_json": str(input_json_path.resolve()),
          "status": "ready",
          "reason": "",
        }
      )

    manifest = {
      "total_samples": len(rows),
      "ready_samples": ready_count,
      "skipped_samples": len(rows) - ready_count,
      "skip_reasons": skip_reasons,
      "records": manifest_records,
    }
    manifest_path = self.out_dir / "_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[OK] boltzgen preprocess 输出: {self.out_dir}")
    print(f"[INFO] boltzgen manifest: {manifest_path}")


def main() -> int:
  parser = argparse.ArgumentParser(description="boltzgen 预处理")
  parser.add_argument("--af3-input", required=True, type=Path)
  parser.add_argument("--dataset-csv", required=True, type=Path)
  parser.add_argument("--out-dir", required=True, type=Path)
  args = parser.parse_args()
  PreProcess(args.af3_input, args.dataset_csv, args.out_dir).preprocess()
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
