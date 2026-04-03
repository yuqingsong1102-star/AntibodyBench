#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Iterable

from evaluation.schema import EVALUATION_LONG_FIELDS


def _read_fasta_sequence(path: Path) -> str:
  if not path.exists():
    return ""
  seq_parts: list[str] = []
  with path.open("r", encoding="utf-8") as f:
    for line in f:
      s = line.strip()
      if not s or s.startswith(">"):
        continue
      seq_parts.append(s)
  return "".join(seq_parts).strip()


def _load_index_by_sample(index_csv: Path) -> dict[str, dict[str, str]]:
  by_sample: dict[str, dict[str, str]] = {}
  if not index_csv.exists():
    return by_sample
  with index_csv.open("r", encoding="utf-8", newline="") as f:
    for row in csv.DictReader(f):
      sid = (row.get("sample_id") or "").strip()
      if sid:
        by_sample[sid] = row
  return by_sample


def _iter_manifest_csvs(native_roots: Iterable[Path]) -> Iterable[tuple[str, str, Path]]:
  for native_root in native_roots:
    run_tag = native_root.name
    if not native_root.exists():
      continue
    for model_dir in sorted(p for p in native_root.iterdir() if p.is_dir()):
      manifest_csv = model_dir / "manifest.csv"
      if manifest_csv.exists():
        yield run_tag, model_dir.name, manifest_csv


def _safe_read_meta(meta_path: Path) -> dict:
  if not meta_path.exists():
    return {}
  try:
    return json.loads(meta_path.read_text(encoding="utf-8"))
  except Exception:
    return {}


def _build_row(
  *,
  run_tag: str,
  model: str,
  manifest_row_idx: int,
  manifest_row: dict[str, str],
  index_row: dict[str, str],
) -> dict[str, str]:
  sample_id = (manifest_row.get("sample_id") or "").strip()
  status = (manifest_row.get("status") or "").strip()
  error_summary = (manifest_row.get("error_summary") or "").strip()
  duration_sec = (manifest_row.get("duration_sec") or "").strip()
  meta_path_str = (manifest_row.get("meta_path") or "").strip()
  sequence_path_str = (manifest_row.get("sequence_path") or "").strip()
  structure_path_str = (manifest_row.get("structure_path") or "").strip()
  meta_path = Path(meta_path_str) if meta_path_str else None
  sequence_path = Path(sequence_path_str) if sequence_path_str else None
  structure_path = Path(structure_path_str) if structure_path_str else None

  meta = _safe_read_meta(meta_path) if meta_path else {}

  pred_sequence = _read_fasta_sequence(sequence_path) if sequence_path else ""
  pred_seq_len = str(len(pred_sequence)) if pred_sequence else ""
  has_sequence = "1" if pred_sequence else "0"
  has_structure = "1" if (structure_path and structure_path.exists() and structure_path.is_file()) else "0"

  out = {k: "" for k in EVALUATION_LONG_FIELDS}
  out["run_id"] = f"{run_tag}:{model}:{sample_id}:{manifest_row_idx}"
  out["run_tag"] = run_tag
  out["model"] = model
  out["manifest_row_idx"] = str(manifest_row_idx)
  out["sample_id"] = sample_id
  out["status"] = status or str(meta.get("status", "")).strip()
  out["error_summary"] = error_summary or str(meta.get("error_summary", "")).strip()
  out["duration_sec"] = duration_sec
  out["meta_path"] = str(meta_path) if meta_path else ""
  out["sequence_path"] = str(sequence_path) if sequence_path else ""
  out["structure_path"] = str(structure_path) if structure_path else ""
  out["pred_sequence"] = pred_sequence
  out["pred_seq_len"] = pred_seq_len
  out["has_structure"] = has_structure
  out["has_sequence"] = has_sequence

  out["reference_complex_path"] = (index_row.get("reference_complex_path") or "").strip()
  out["reference_complex_status"] = (index_row.get("reference_complex_status") or "").strip()
  out["antigen_chain_id"] = (index_row.get("antigen_chain") or "").strip()
  out["antibody_chain_id"] = (index_row.get("antibody_chain") or "").strip()
  out["cdr_h3_start"] = (index_row.get("cdr_h3_start") or "").strip()
  out["cdr_h3_end"] = (index_row.get("cdr_h3_end") or "").strip()
  out["cdr_h3_status"] = (index_row.get("cdr_h3_status") or "").strip()

  # 该行先只负责入库，指标阶段再计算这些字段。
  out["aa_valid_ratio"] = ""
  out["cys_count"] = ""
  out["motif_h3_len_in_8_20"] = ""
  out["cdr_h3_rmsd"] = ""
  out["cdr_h3_atom_count"] = ""
  out["struct_eval_eligible"] = ""
  out["metric_error"] = ""
  return out


def build_evaluation_long(native_roots: list[Path], index_csv: Path, out_csv: Path) -> Path:
  index_by_sample = _load_index_by_sample(index_csv)
  rows: list[dict[str, str]] = []

  for run_tag, model, manifest_csv in _iter_manifest_csvs(native_roots):
    with manifest_csv.open("r", encoding="utf-8", newline="") as f:
      reader = csv.DictReader(f)
      for i, mrow in enumerate(reader, start=1):
        sid = (mrow.get("sample_id") or "").strip()
        idx_row = index_by_sample.get(sid, {})
        rows.append(
          _build_row(
            run_tag=run_tag,
            model=model,
            manifest_row_idx=i,
            manifest_row=mrow,
            index_row=idx_row,
          )
        )

  out_csv.parent.mkdir(parents=True, exist_ok=True)
  with out_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=EVALUATION_LONG_FIELDS)
    writer.writeheader()
    writer.writerows(rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="将 native_predictions 输出拼接为统一 evaluation_long.csv")
  parser.add_argument(
    "--native-root",
    action="append",
    type=Path,
    help="native 输出目录（可重复指定）。默认自动读取 real+smoke。",
  )
  parser.add_argument(
    "--index-csv",
    type=Path,
    default=Path("inputs/dataset_index_h3_annotated.csv"),
    help="样本索引 CSV（建议使用包含 cdr_h3_status 的版本）",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/evaluation_long.csv"),
    help="输出长表路径",
  )
  args = parser.parse_args()

  native_roots = args.native_root or [
    Path("outputs/native_predictions_real"),
    Path("outputs/native_predictions_smoke"),
  ]
  out_csv = build_evaluation_long(native_roots=native_roots, index_csv=args.index_csv, out_csv=args.out_csv)
  print(f"[OK] 写出评估长表: {out_csv}")
  print(f"[INFO] native_roots={','.join(str(p) for p in native_roots)}")
  print(f"[INFO] index_csv={args.index_csv}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

