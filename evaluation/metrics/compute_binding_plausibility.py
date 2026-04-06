#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import csv
import math
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser
from Bio.PDB.SASA import ShrakeRupley

from evaluation.schema import EXTERNAL_SEED_FIELDS


AA20 = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC_RES3 = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "PRO"}

CONTACT_CUTOFF = 5.0
CLASH_CUTOFF = 2.0
SURFACE_SASA_THRESHOLD = 15.0

# 9.1 初始占位阈值
MIN_CONTACT_PAIRS = 15
MIN_INTERFACE_BSA = 400.0
MAX_CLASH_RATIO = 0.05
MIN_HOTSPOT_CONTACTS = 1
MIN_PDOCKQ2 = 0.23
MIN_EXTERNAL_IPTM = 0.60
MAX_EXTERNAL_IPAE = 10.0


def _load_csv(path: Path) -> list[dict[str, str]]:
  if not path.exists():
    return []
  with path.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=EXTERNAL_SEED_FIELDS)
    writer.writeheader()
    writer.writerows(rows)


def _round_or_blank(value: float | None) -> str:
  if value is None:
    return ""
  return str(round(value, 6))


def _to_bool_str(flag: bool) -> str:
  return "1" if flag else "0"


def _to_float(text: str) -> float | None:
  value = (text or "").strip()
  if not value:
    return None
  try:
    return float(value)
  except Exception:
    return None


def _parse_structure(path: Path):
  if path.suffix.lower() in {".cif", ".mmcif"}:
    parser = MMCIFParser(QUIET=True)
  else:
    parser = PDBParser(QUIET=True)
  return parser.get_structure("s", str(path))


def _first_model(structure):
  return next(iter(structure))


def _get_chain(model, chain_id: str):
  if not chain_id:
    return None
  for chain in model.get_chains():
    if chain.id == chain_id:
      return chain
  return None


def _res_key(residue) -> tuple[int, str]:
  het, resseq, icode = residue.get_id()
  return int(resseq), str(icode or " ")


def _res_label(chain_id: str, residue) -> str:
  resseq, icode = _res_key(residue)
  icode_text = "" if icode == " " else icode
  return f"{chain_id}{resseq}{icode_text}"


def _iter_heavy_atoms(chain):
  atoms = []
  for atom in chain.get_atoms():
    element = (getattr(atom, "element", "") or "").strip().upper()
    if element == "H":
      continue
    atoms.append(atom)
  return atoms


def _extract_single_chain_structure(structure, chain_id: str):
  structure_copy = copy.deepcopy(structure)
  model = _first_model(structure_copy)
  for chain in list(model.get_chains()):
    if chain.id != chain_id:
      model.detach_child(chain.id)
  return structure_copy


def _sum_chain_atom_sasa(chain) -> float:
  return float(sum(float(getattr(atom, "sasa", 0.0)) for atom in _iter_heavy_atoms(chain)))


def _compute_interface_bsa(structure, target_chain_id: str, antibody_chain_id: str) -> tuple[float | None, float | None]:
  try:
    sr = ShrakeRupley()
    full_structure = copy.deepcopy(structure)
    sr.compute(full_structure, level="A")
    full_model = _first_model(full_structure)
    target_bound = _sum_chain_atom_sasa(_get_chain(full_model, target_chain_id))
    antibody_bound = _sum_chain_atom_sasa(_get_chain(full_model, antibody_chain_id))

    target_only = _extract_single_chain_structure(structure, target_chain_id)
    antibody_only = _extract_single_chain_structure(structure, antibody_chain_id)
    sr.compute(target_only, level="A")
    sr.compute(antibody_only, level="A")
    target_unbound = _sum_chain_atom_sasa(_get_chain(_first_model(target_only), target_chain_id))
    antibody_unbound = _sum_chain_atom_sasa(_get_chain(_first_model(antibody_only), antibody_chain_id))

    target_bsa = max(0.0, target_unbound - target_bound)
    antibody_bsa = max(0.0, antibody_unbound - antibody_bound)
    interface_bsa = 0.5 * (target_bsa + antibody_bsa)
    binder_fraction = None
    if antibody_unbound > 0:
      binder_fraction = (antibody_bsa / antibody_unbound) * 100.0
    return interface_bsa, binder_fraction
  except Exception:
    return None, None


def _compute_surface_hydrophobicity(structure, antibody_chain_id: str) -> float | None:
  try:
    sr = ShrakeRupley()
    antibody_only = _extract_single_chain_structure(structure, antibody_chain_id)
    sr.compute(antibody_only, level="R")
    chain = _get_chain(_first_model(antibody_only), antibody_chain_id)
    surface_total = 0
    surface_hydrophobic = 0
    for residue in chain.get_residues():
      het, _resseq, _icode = residue.get_id()
      if het != " ":
        continue
      sasa = float(getattr(residue, "sasa", 0.0))
      if sasa < SURFACE_SASA_THRESHOLD:
        continue
      surface_total += 1
      if residue.get_resname().strip().upper() in HYDROPHOBIC_RES3:
        surface_hydrophobic += 1
    if surface_total == 0:
      return None
    return surface_hydrophobic / surface_total
  except Exception:
    return None


def _compute_pdockq2(target_chain, antibody_chain) -> tuple[float | None, float | None]:
  target_coords = []
  target_plddt = []
  antibody_coords = []
  antibody_plddt = []
  for chain, coords, plddt in [
    (target_chain, target_coords, target_plddt),
    (antibody_chain, antibody_coords, antibody_plddt),
  ]:
    for residue in chain.get_residues():
      het, _resseq, _icode = residue.get_id()
      if het != " ":
        continue
      atom = residue["CB"] if "CB" in residue else (residue["CA"] if "CA" in residue else None)
      if atom is None:
        continue
      coords.append(np.asarray(atom.get_coord(), dtype=float))
      plddt.append(float(atom.get_bfactor()))
  if not target_coords or not antibody_coords:
    return None, None
  target_coords_arr = np.asarray(target_coords, dtype=float)
  antibody_coords_arr = np.asarray(antibody_coords, dtype=float)
  dists = np.sqrt(np.sum((target_coords_arr[:, None, :] - antibody_coords_arr[None, :, :]) ** 2, axis=-1))
  contact_mask = dists <= 8.0
  n_contacts = int(np.count_nonzero(contact_mask))
  if n_contacts == 0:
    return 0.0, 0.0
  target_interface = np.where(contact_mask.any(axis=1))[0]
  antibody_interface = np.where(contact_mask.any(axis=0))[0]
  plddt_values = [target_plddt[i] for i in target_interface] + [antibody_plddt[i] for i in antibody_interface]
  if not plddt_values or max(abs(v) for v in plddt_values) < 1e-6:
    return None, None
  avg_interface_plddt = float(np.mean(plddt_values))
  x = avg_interface_plddt * math.log10(n_contacts)
  # 这里沿用 pDockQ 逻辑作为 pDockQ2 占位估计。
  pdockq2 = 0.724 / (1.0 + math.exp(-0.052 * (x - 152.611))) + 0.018
  return pdockq2, avg_interface_plddt


def _mean_binder_plddt(chain) -> float | None:
  values = []
  for residue in chain.get_residues():
    het, _resseq, _icode = residue.get_id()
    if het != " ":
      continue
    atom = residue["CA"] if "CA" in residue else None
    if atom is None:
      continue
    values.append(float(atom.get_bfactor()))
  if not values or max(abs(v) for v in values) < 1e-6:
    return None
  return float(np.mean(values))


def _compute_interface_metrics(structure, target_chain_id: str, antibody_chain_id: str, hotspot_labels: set[str]):
  model = _first_model(structure)
  target_chain = _get_chain(model, target_chain_id)
  antibody_chain = _get_chain(model, antibody_chain_id)
  if target_chain is None or antibody_chain is None:
    raise KeyError("missing_chain")

  target_atoms = _iter_heavy_atoms(target_chain)
  antibody_atoms = _iter_heavy_atoms(antibody_chain)
  atom_search = NeighborSearch(target_atoms + antibody_atoms)
  target_atom_ids = {id(atom) for atom in target_atoms}

  contact_res_pairs: set[tuple[str, str]] = set()
  target_contact_residues: set[str] = set()
  antibody_contact_residues: set[str] = set()
  contact_atom_pairs = 0
  clash_atom_pairs = 0

  for antibody_atom in antibody_atoms:
    nearby_atoms = atom_search.search(antibody_atom.get_coord(), CONTACT_CUTOFF, level="A")
    for target_atom in nearby_atoms:
      if id(target_atom) not in target_atom_ids:
        continue
      target_residue = target_atom.get_parent()
      antibody_residue = antibody_atom.get_parent()
      target_label = _res_label(target_chain_id, target_residue)
      antibody_label = _res_label(antibody_chain_id, antibody_residue)
      contact_res_pairs.add((target_label, antibody_label))
      target_contact_residues.add(target_label)
      antibody_contact_residues.add(antibody_label)
      contact_atom_pairs += 1
      if float(target_atom - antibody_atom) < CLASH_CUTOFF:
        clash_atom_pairs += 1

  interface_residue_count = len(antibody_contact_residues)
  antigen_contact_residue_count = len(target_contact_residues)
  interface_contact_pair_count = len(contact_res_pairs)
  hotspot_contact_count = len(target_contact_residues & hotspot_labels)
  interface_clash_ratio = (clash_atom_pairs / contact_atom_pairs) if contact_atom_pairs else None

  antibody_resseqs = sorted({int(_res_key(res)[0]) for res in (atom.get_parent() for atom in antibody_atoms) if res.get_id()[0] == " "})
  segment_count = 0
  prev = None
  for resseq in antibody_resseqs:
    if prev is None or resseq > prev + 1:
      segment_count += 1
    prev = resseq
  fragmentation_index = (segment_count / interface_residue_count) if interface_residue_count else None

  interface_bsa, binder_fraction = _compute_interface_bsa(structure, target_chain_id, antibody_chain_id)
  pdockq2, avg_interface_plddt = _compute_pdockq2(target_chain, antibody_chain)

  return {
    "interface_residue_count": float(interface_residue_count),
    "antigen_contact_residue_count": float(antigen_contact_residue_count),
    "interface_contact_pair_count": float(interface_contact_pair_count),
    "interface_bsa": interface_bsa,
    "interface_binder_fraction": binder_fraction,
    "interface_clash_count": float(clash_atom_pairs),
    "interface_clash_ratio": interface_clash_ratio,
    "hotspot_contact_count": float(hotspot_contact_count),
    "interface_fragmentation_index": fragmentation_index,
    "external_plddt_binder": _mean_binder_plddt(antibody_chain),
    "pdockq2": pdockq2,
    "avg_interface_plddt": avg_interface_plddt,
    "surface_hydrophobicity": _compute_surface_hydrophobicity(structure, antibody_chain_id),
    # Phase 1 占位字段
    "interface_shape_complementarity": None,
    "interface_packstat": None,
  }


def _score_row(row: dict[str, str]) -> dict[str, str]:
  out = {k: (row.get(k, "") or "") for k in EXTERNAL_SEED_FIELDS}
  seq = (out.get("pred_sequence") or "").strip().upper()
  valid_ratio = (sum(1 for aa in seq if aa in AA20) / len(seq)) if seq else None
  out["sequence_valid_ratio"] = _round_or_blank(valid_ratio)
  out["sequence_length"] = out.get("sequence_length") or (str(len(seq)) if seq else "")
  out["pred_seq_len"] = out.get("pred_seq_len") or out["sequence_length"]
  out["has_valid_sequence"] = _to_bool_str(bool(seq) and valid_ratio == 1.0)

  structure_path = Path((out.get("external_structure_path") or "").strip()) if (out.get("external_structure_path") or "").strip() else None
  has_external_structure = bool(structure_path and structure_path.exists())
  out["has_external_structure"] = _to_bool_str(has_external_structure)
  out["external_prediction_ok"] = _to_bool_str((out.get("external_status") or "").strip() == "ok" and has_external_structure)

  hotspot_labels = {part.strip() for part in (out.get("hotspot_string") or "").split(",") if part.strip()}
  target_chain_id = (out.get("target_chain_id") or "").strip()
  antibody_chain_id = (out.get("designed_antibody_chain_id") or "").strip()

  chain_parse_ok = False
  severe_clash_flag = False
  metric_error_code = ""
  metric_error = ""

  if has_external_structure and out["external_prediction_ok"] == "1":
    try:
      structure = _parse_structure(structure_path)
      model = _first_model(structure)
      if _get_chain(model, target_chain_id) is None or _get_chain(model, antibody_chain_id) is None:
        raise KeyError("missing_chain")
      chain_parse_ok = True
      metrics = _compute_interface_metrics(structure, target_chain_id, antibody_chain_id, hotspot_labels)
      for key, value in metrics.items():
        out[key] = _round_or_blank(value)
      clash_ratio = metrics.get("interface_clash_ratio")
      clash_count = metrics.get("interface_clash_count")
      severe_clash_flag = ((clash_ratio is not None and clash_ratio > 0.1) or (clash_count is not None and clash_count > 20))
      out["surface_hydrophobicity"] = _round_or_blank(metrics.get("surface_hydrophobicity"))
      out["pdockq2"] = _round_or_blank(metrics.get("pdockq2"))
      out["external_plddt_binder"] = _round_or_blank(metrics.get("external_plddt_binder"))
    except KeyError as exc:
      metric_error_code = str(exc).strip("'\"")
      metric_error = metric_error_code
    except Exception as exc:
      metric_error_code = f"binding_metric_failed:{type(exc).__name__}"
      metric_error = metric_error_code

  out["chain_parse_ok"] = _to_bool_str(chain_parse_ok)
  out["severe_clash_flag"] = _to_bool_str(severe_clash_flag)
  viability_pass = (
    out["has_valid_sequence"] == "1"
    and out["has_external_structure"] == "1"
    and out["external_prediction_ok"] == "1"
    and out["chain_parse_ok"] == "1"
  )
  out["viability_pass"] = _to_bool_str(viability_pass)
  if not viability_pass:
    reasons = []
    if out["has_valid_sequence"] != "1":
      reasons.append("invalid_sequence")
    if out["has_external_structure"] != "1":
      reasons.append("missing_external_structure")
    if out["external_prediction_ok"] != "1":
      reasons.append("external_prediction_failed")
    if out["chain_parse_ok"] != "1":
      reasons.append("chain_parse_failed")
    out["viability_fail_reason"] = ",".join(reasons) if reasons else "unknown"
  else:
    out["viability_fail_reason"] = ""

  binding_checks = []
  interface_contact_pairs = _to_float(out.get("interface_contact_pair_count", ""))
  interface_bsa = _to_float(out.get("interface_bsa", ""))
  interface_clash_ratio = _to_float(out.get("interface_clash_ratio", ""))
  hotspot_contact_count = _to_float(out.get("hotspot_contact_count", ""))
  pdockq2 = _to_float(out.get("pdockq2", ""))
  external_iptm = _to_float(out.get("external_iptm", ""))
  external_ipae = _to_float(out.get("external_ipae", ""))

  if interface_contact_pairs is not None:
    binding_checks.append(interface_contact_pairs >= MIN_CONTACT_PAIRS)
  if interface_bsa is not None:
    binding_checks.append(interface_bsa >= MIN_INTERFACE_BSA)
  if interface_clash_ratio is not None:
    binding_checks.append(interface_clash_ratio <= MAX_CLASH_RATIO)
  if hotspot_labels and hotspot_contact_count is not None:
    binding_checks.append(hotspot_contact_count >= MIN_HOTSPOT_CONTACTS)
  if pdockq2 is not None:
    binding_checks.append(pdockq2 > MIN_PDOCKQ2)
  if external_iptm is not None:
    binding_checks.append(external_iptm >= MIN_EXTERNAL_IPTM)
  if external_ipae is not None:
    binding_checks.append(external_ipae <= MAX_EXTERNAL_IPAE)

  binding_pass = viability_pass and bool(binding_checks) and all(binding_checks)
  out["binding_plausibility_pass"] = _to_bool_str(binding_pass)
  out["strong_binder_proxy_pass"] = _to_bool_str(binding_pass and (pdockq2 or 0.0) >= 0.49 and (interface_bsa or 0.0) >= 600.0)

  out["metric_error_code"] = metric_error_code
  out["metric_error"] = metric_error
  return out


def compute_binding_plausibility(in_csv: Path, out_csv: Path) -> Path:
  rows = _load_csv(in_csv)
  out_rows = [_score_row(row) for row in rows]
  _write_csv(out_csv, out_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="计算 Viability 与 Binding Plausibility 评分层")
  parser.add_argument(
    "--in-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="seed 级主表",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="输出 seed 级主表（可覆盖输入）",
  )
  args = parser.parse_args()
  out = compute_binding_plausibility(args.in_csv, args.out_csv)
  print(f"[OK] binding plausibility: {out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
