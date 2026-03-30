#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any


class PostProcess:
  def __init__(self, dataset_csv: Path, prediction_dir: Path, out_dir: Path) -> None:
    self.dataset_csv = dataset_csv
    self.prediction_dir = prediction_dir
    self.out_dir = out_dir

  def _parse_optional_int(self, row: dict[str, str], key: str) -> int | None:
    v = (row.get(key) or "").strip()
    if not v:
      return None
    try:
      return int(v)
    except ValueError:
      return None

  def _default_ref_path(self, pdb_id_raw: str) -> Path | None:
    pdb_base = (pdb_id_raw or "").split("-")[0].strip()
    if not pdb_base:
      return None
    root_dir = Path(__file__).resolve().parents[2]
    cand = root_dir / "data" / "raw" / "raw_pdbs" / f"{pdb_base}.pdb"
    return cand if cand.exists() else None

  def postprocess(self) -> None:
    self.out_dir.mkdir(parents=True, exist_ok=True)
    cdr_dir = self.out_dir / "cdr_regions"
    cdr_dir.mkdir(parents=True, exist_ok=True)

    root_dir = Path(__file__).resolve().parents[2]
    if str(root_dir) not in sys.path:
      sys.path.insert(0, str(root_dir))
    from evaluation_tools.cdr_rmsd import compute_antigen_aligned_cdr_rmsd, extract_cdr_backbone_pdb

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

      pdb_id = str(row.get("pdb_id", "")).strip()
      ref_path_raw = str(row.get("reference_structure_path", "")).strip()
      ref_path = Path(ref_path_raw) if ref_path_raw else self._default_ref_path(pdb_id)
      if ref_path is None or not ref_path.exists():
        continue

      antigen_chain = (row.get("antigen_chain") or "").strip() or "A"
      heavy_chain = (row.get("antibody_heavy_chain") or row.get("antibody_chain") or "").strip() or ""
      light_chain = (row.get("antibody_light_chain") or row.get("antibody_chain") or "").strip() or ""

      h3_start = self._parse_optional_int(row, "cdr_h3_start")
      h3_end = self._parse_optional_int(row, "cdr_h3_end")
      l3_start = self._parse_optional_int(row, "cdr_l3_start")
      l3_end = self._parse_optional_int(row, "cdr_l3_end")

      cdr_h3_rmsd: str = ""
      cdr_l3_rmsd: str = ""

      if h3_start is not None and h3_end is not None and heavy_chain:
        res = compute_antigen_aligned_cdr_rmsd(
          pred_path,
          ref_path,
          antigen_chain_id=antigen_chain,
          antibody_chain_id=heavy_chain,
          cdr_start=h3_start,
          cdr_end=h3_end,
        )
        if res is not None:
          cdr_h3_rmsd = f"{res.rmsd:.6f}"
        try:
          extract_cdr_backbone_pdb(
            pred_path,
            cdr_dir / f"{sample_id}_cdr_h3.pdb",
            chain_id=heavy_chain,
            start=h3_start,
            end=h3_end,
          )
        except Exception:
          pass

      if l3_start is not None and l3_end is not None and light_chain:
        res = compute_antigen_aligned_cdr_rmsd(
          pred_path,
          ref_path,
          antigen_chain_id=antigen_chain,
          antibody_chain_id=light_chain,
          cdr_start=l3_start,
          cdr_end=l3_end,
        )
        if res is not None:
          cdr_l3_rmsd = f"{res.rmsd:.6f}"
        try:
          extract_cdr_backbone_pdb(
            pred_path,
            cdr_dir / f"{sample_id}_cdr_l3.pdb",
            chain_id=light_chain,
            start=l3_start,
            end=l3_end,
          )
        except Exception:
          pass

      prediction_path = str(pred_path.resolve())

      rows.append(
        {
          "pdb_id": pdb_id,
          "seed": 0,
          "sample": sample_id,
          "ranking_score": "",
          "prediction_path": prediction_path,
          "cdr_h3_rmsd": cdr_h3_rmsd,
          "cdr_l3_rmsd": cdr_l3_rmsd,
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
    print(f"[OK] BindCraft postprocess 输出: {out_csv}")


def main() -> int:
  parser = argparse.ArgumentParser(description="BindCraft 后处理")
  parser.add_argument("--dataset-csv", required=True, type=Path)
  parser.add_argument("--prediction-dir", required=True, type=Path)
  parser.add_argument("--out-dir", required=True, type=Path)
  args = parser.parse_args()
  PostProcess(args.dataset_csv, args.prediction_dir, args.out_dir).postprocess()
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
