#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any


AA3_TO_1 = {
  "ALA": "A",
  "ARG": "R",
  "ASN": "N",
  "ASP": "D",
  "CYS": "C",
  "GLN": "Q",
  "GLU": "E",
  "GLY": "G",
  "HIS": "H",
  "ILE": "I",
  "LEU": "L",
  "LYS": "K",
  "MET": "M",
  "PHE": "F",
  "PRO": "P",
  "SER": "S",
  "THR": "T",
  "TRP": "W",
  "TYR": "Y",
  "VAL": "V",
  "SEC": "U",
  "PYL": "O",
}


@dataclass
class PickResult:
  structure_path: Path | None
  sequence: str
  selection_rule: str
  sequence_source: str
  notes: list[str]


def _first_existing(paths: list[Path]) -> Path | None:
  for p in paths:
    if p.exists():
      return p
  return None


def _first_glob(root: Path, patterns: list[str]) -> Path | None:
  for pattern in patterns:
    matches = sorted(root.glob(pattern))
    if matches:
      return matches[0]
  return None


def _read_csv_first_row(path: Path) -> dict[str, str]:
  if not path.exists():
    return {}
  with path.open("r", encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
      return {str(k): str(v) for k, v in row.items()}
  return {}


def _sequence_from_row(row: dict[str, str]) -> tuple[str, str]:
  for key, value in row.items():
    key_norm = key.strip().lower()
    if "sequence" in key_norm or key_norm in {"seq", "aa"}:
      seq = (value or "").strip().replace(" ", "")
      if seq:
        return seq, f"csv:{key}"
  return "", ""


def _infer_binder_chain(sample_id: str) -> tuple[str | None, str | None]:
  parts = sample_id.split("_")
  if len(parts) >= 3:
    antibody_chain = parts[1].strip() or None
    antigen_chain = parts[2].strip() or None
    return antibody_chain, antigen_chain
  return None, None


def _sequence_from_pdb(structure_path: Path, preferred_chain: str | None, avoid_chain: str | None) -> tuple[str, str]:
  chain_residues: dict[str, list[str]] = {}
  seen: set[tuple[str, str, str, str]] = set()
  with structure_path.open("r", encoding="utf-8", errors="ignore") as f:
    for line in f:
      if not (line.startswith("ATOM") or line.startswith("HETATM")):
        continue
      if len(line) < 27:
        continue
      chain_id = line[21].strip() or "_"
      resn = line[17:20].strip().upper()
      resseq = line[22:26].strip()
      icode = line[26].strip() if len(line) > 26 else ""
      atom_name = line[12:16].strip()
      # 用每个残基第一个原子去重，避免重复计数
      key = (chain_id, resseq, icode, atom_name[:1])
      residue_key = (chain_id, resseq, icode, "R")
      if residue_key in seen:
        continue
      if atom_name not in {"CA", "C1'", "C1*"}:
        # 没有 CA 时，依旧允许首个原子代表该残基
        if key in seen:
          continue
      seen.add(residue_key)
      aa = AA3_TO_1.get(resn, "X")
      chain_residues.setdefault(chain_id, []).append(aa)

  if not chain_residues:
    return "", ""

  if preferred_chain and preferred_chain in chain_residues:
    return "".join(chain_residues[preferred_chain]), f"pdb_chain:{preferred_chain}"

  if avoid_chain:
    candidates = {k: v for k, v in chain_residues.items() if k != avoid_chain}
    if candidates:
      best_chain = max(candidates.keys(), key=lambda c: len(candidates[c]))
      return "".join(candidates[best_chain]), f"pdb_chain:{best_chain}"

  best_chain = max(chain_residues.keys(), key=lambda c: len(chain_residues[c]))
  return "".join(chain_residues[best_chain]), f"pdb_chain:{best_chain}"


def _sequence_from_structure(path: Path, sample_id: str) -> tuple[str, str]:
  if path.suffix.lower() != ".pdb":
    return "", ""
  preferred_chain, avoid_chain = _infer_binder_chain(sample_id)
  return _sequence_from_pdb(path, preferred_chain=preferred_chain, avoid_chain=avoid_chain)


def _resolve_candidate(path_str: str, base_dir: Path, native_root: Path) -> Path | None:
  raw = (path_str or "").strip()
  if not raw:
    return None
  p = Path(raw)
  if p.is_absolute() and p.exists():
    return p
  local = (base_dir / raw).resolve()
  if local.exists():
    return local
  alt = (native_root / raw).resolve()
  if alt.exists():
    return alt
  return None


def pick_rfantibody(native_root: Path, sample_id: str) -> PickResult:
  structure = _first_glob(
    native_root,
    [
      "**/qv/**/*.pdb",
      "**/qv/**/*.cif",
      "**/*qv*/*.pdb",
      "**/*qv*/*.cif",
      "**/*.pdb",
      "**/*.cif",
    ],
  )
  notes: list[str] = []
  seq = ""
  source = ""
  if structure:
    seq, source = _sequence_from_structure(structure, sample_id)
    if not seq:
      notes.append("sequence_not_found_in_rfantibody_outputs")
  return PickResult(structure, seq, "rfantibody_qv_first", source, notes)


def pick_germinal(native_root: Path, sample_id: str) -> PickResult:
  notes: list[str] = []
  csv_path = _first_glob(native_root, ["**/accepted/designs.csv", "**/designs.csv"])
  row = _read_csv_first_row(csv_path) if csv_path else {}
  structure: Path | None = None
  if row and csv_path:
    for key in row:
      key_norm = key.lower()
      if "path" in key_norm or "pdb" in key_norm or "structure" in key_norm:
        structure = _resolve_candidate(row[key], csv_path.parent, native_root)
        if structure:
          break
  if not structure:
    structure = _first_glob(native_root, ["**/accepted/structures/*.pdb", "**/accepted/structures/*.cif"])

  seq, source = _sequence_from_row(row)
  if not seq and structure:
    seq, source = _sequence_from_structure(structure, sample_id)
    if seq:
      notes.append("sequence_fallback_from_structure")
  if not seq:
    notes.append("sequence_not_found_in_germinal_outputs")
  return PickResult(structure, seq, "germinal_designs_csv_first_row", source, notes)


def pick_bindcraft(native_root: Path, sample_id: str) -> PickResult:
  notes: list[str] = []
  csv_path = _first_glob(native_root, ["**/final_design_stats.csv", "**/*design_stats*.csv"])
  row = _read_csv_first_row(csv_path) if csv_path else {}
  structure: Path | None = None
  if row and csv_path:
    for key in row:
      key_norm = key.lower()
      if "path" in key_norm or "pdb" in key_norm or "structure" in key_norm:
        structure = _resolve_candidate(row[key], csv_path.parent, native_root)
        if structure:
          break
  if not structure:
    structure = _first_glob(native_root, ["**/Accepted/*.pdb", "**/Accepted/*.cif", "**/accepted/*.pdb"])

  seq, source = _sequence_from_row(row)
  if not seq and structure:
    seq, source = _sequence_from_structure(structure, sample_id)
    if seq:
      notes.append("sequence_fallback_from_structure")
  if not seq:
    notes.append("sequence_not_found_in_bindcraft_outputs")
  return PickResult(structure, seq, "bindcraft_final_stats_first_row", source, notes)


def _rank1_from_dir(root: Path) -> Path | None:
  ranked = sorted(root.glob("**/*rank*1*.cif")) + sorted(root.glob("**/*rank*1*.pdb"))
  if ranked:
    return ranked[0]
  any_struct = sorted(root.glob("**/*.cif")) + sorted(root.glob("**/*.pdb"))
  if any_struct:
    return any_struct[0]
  return None


def pick_boltzgen(native_root: Path, sample_id: str) -> PickResult:
  notes: list[str] = []
  final_dir = _first_glob(native_root, ["**/final_ranked_designs/final_*_designs", "**/final_ranked_designs"])
  inter_dir = _first_glob(
    native_root,
    ["**/intermediate_ranked_designs/intermediate_*_designs", "**/intermediate_ranked_designs"],
  )
  structure = None
  rule = "boltzgen_final_rank1_then_intermediate_rank1"
  if final_dir:
    structure = _rank1_from_dir(final_dir)
  if not structure and inter_dir:
    structure = _rank1_from_dir(inter_dir)

  seq = ""
  source = ""
  metrics_csv = _first_glob(native_root, ["**/*metrics*.csv", "**/*summary*.csv"])
  if metrics_csv:
    row = _read_csv_first_row(metrics_csv)
    seq, source = _sequence_from_row(row)
  if not seq and structure:
    seq, source = _sequence_from_structure(structure, sample_id)
    if seq:
      notes.append("sequence_fallback_from_structure")
  if not seq:
    notes.append("sequence_not_found_in_boltzgen_outputs")

  return PickResult(structure, seq, rule, source, notes)


def write_fasta(path: Path, header: str, sequence: str) -> None:
  path.write_text(f">{header}\n{sequence}\n", encoding="utf-8")


def main() -> int:
  parser = argparse.ArgumentParser(description="Collect Top-1 outputs into unified contract.")
  parser.add_argument("--model", required=True, type=str)
  parser.add_argument("--sample-id", required=True, type=str)
  parser.add_argument("--sample-input-dir", required=True, type=Path)
  parser.add_argument("--sample-output-dir", required=True, type=Path)
  parser.add_argument("--runner-exit-code", required=True, type=int)
  args = parser.parse_args()

  model = args.model
  sample_id = args.sample_id
  sample_out = args.sample_output_dir
  sample_out.mkdir(parents=True, exist_ok=True)

  native_root = sample_out / "native_run"
  top1_meta = sample_out / "top1_meta.json"
  status = "ok"
  error_summary = ""
  notes: list[str] = []

  picker_map: dict[str, Any] = {
    "RFantibody": pick_rfantibody,
    "germinal": pick_germinal,
    "BindCraft": pick_bindcraft,
    "boltzgen": pick_boltzgen,
  }
  picker = picker_map.get(model)
  if picker is None:
    status = "failed"
    error_summary = f"unsupported_model:{model}"
    result = PickResult(None, "", "unsupported", "", [])
  elif not native_root.exists():
    status = "failed"
    error_summary = "native_run_dir_not_found"
    result = PickResult(None, "", "native_output_missing", "", [])
  else:
    result = picker(native_root, sample_id)
    notes.extend(result.notes)

  top1_structure_path = ""
  top1_sequence_path = ""

  if result.structure_path and result.structure_path.exists():
    ext = result.structure_path.suffix.lower() or ".pdb"
    dst_structure = sample_out / f"top1_structure{ext}"
    shutil.copy2(result.structure_path, dst_structure)
    top1_structure_path = str(dst_structure)
  else:
    status = "failed"
    if not error_summary:
      error_summary = "top1_structure_not_found"

  sequence = (result.sequence or "").strip()
  if sequence:
    dst_fasta = sample_out / "top1_sequence.fasta"
    write_fasta(dst_fasta, f"{model}|{sample_id}|top1", sequence)
    top1_sequence_path = str(dst_fasta)
  else:
    if status != "failed":
      notes.append("empty_top1_sequence")
    dst_fasta = sample_out / "top1_sequence.fasta"
    write_fasta(dst_fasta, f"{model}|{sample_id}|top1", "")
    top1_sequence_path = str(dst_fasta)

  if args.runner_exit_code != 0 and status != "failed":
    notes.append(f"runner_nonzero_exit:{args.runner_exit_code}")
  if args.runner_exit_code != 0 and status == "failed" and not error_summary:
    error_summary = f"runner_exit_nonzero:{args.runner_exit_code}"

  meta = {
    "status": status,
    "model": model,
    "sample_id": sample_id,
    "runner_exit_code": args.runner_exit_code,
    "selection_rule": result.selection_rule,
    "source_structure": str(result.structure_path) if result.structure_path else "",
    "source_sequence": result.sequence_source,
    "top1_structure": top1_structure_path,
    "top1_sequence": top1_sequence_path,
    "error_summary": error_summary,
    "notes": notes,
    "sample_input_dir": str(args.sample_input_dir),
    "sample_output_dir": str(sample_out),
  }
  top1_meta.write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
