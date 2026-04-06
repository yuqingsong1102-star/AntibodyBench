#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import MMCIFParser, NeighborSearch, PDBIO, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.SeqUtils import seq1


BINDER_KINDS = {"heavy", "light", "binder"}


@dataclass
class ChainInfo:
  chain_id: str
  molecule_name: str
  kind: str


@dataclass
class SampleResult:
  sample_id: str
  status: str
  note: str
  antigen_chain_used: str
  antibody_chain_used: str
  hotspot_count: int
  hotspots: list[str]
  target_pdb_repaired: bool


def _parse_structure(path: Path):
  if path.suffix.lower() in {".cif", ".mmcif"}:
    parser = MMCIFParser(QUIET=True)
  else:
    parser = PDBParser(QUIET=True)
  return parser.get_structure("complex", str(path))


def _as_list(value) -> list[str]:
  if value is None:
    return []
  if isinstance(value, list):
    return [str(x) for x in value]
  return [str(value)]


def _normalize_name(text: str) -> str:
  return re.sub(r"\s+", " ", (text or "").upper()).strip()


def _classify_chain_kind(molecule_name: str) -> str:
  name = _normalize_name(molecule_name)
  if not name:
    return "other"
  if "LIGHT CHAIN" in name or "FAB LIGHT" in name or re.search(r"(^|[^A-Z])VL([^A-Z]|$)", name):
    return "light"
  if "HEAVY CHAIN" in name or "FAB HEAVY" in name or re.search(r"(^|[^A-Z])VH([^A-Z]|$)", name):
    return "heavy"
  if "VHH" in name or "NANOBODY" in name or "SCFV" in name or "ANTIBODY" in name or "FAB" in name or "IMMUNOGLOBULIN" in name:
    return "binder"
  return "other"


def _parse_pdb_chain_info(path: Path) -> dict[str, ChainInfo]:
  out: dict[str, ChainInfo] = {}
  if path.suffix.lower() not in {".pdb", ".ent"}:
    return out

  molecules: dict[str, dict[str, str]] = {}
  current_mol_id: str | None = None
  current_key: str | None = None
  with path.open("r", encoding="utf-8", errors="ignore") as f:
    for line in f:
      if not line.startswith("COMPND"):
        continue
      text = line[10:].strip()
      if ":" in text:
        key, value = text.split(":", 1)
        current_key = key.strip().lower()
        value = value.strip().rstrip(";")
        if current_key == "mol_id":
          current_mol_id = value
          molecules.setdefault(current_mol_id, {"molecule": "", "chains": ""})
        elif current_mol_id is not None:
          field = "chains" if current_key == "chain" else current_key
          molecules[current_mol_id][field] = (molecules[current_mol_id].get(field, "") + " " + value).strip()
      elif current_mol_id is not None and current_key in {"molecule", "chain"}:
        field = "chains" if current_key == "chain" else current_key
        molecules[current_mol_id][field] = (molecules[current_mol_id].get(field, "") + " " + text.rstrip(";")).strip()

  for payload in molecules.values():
    molecule_name = _normalize_name(payload.get("molecule", ""))
    kind = _classify_chain_kind(molecule_name)
    chains = [x.strip() for x in re.split(r"[;, ]+", payload.get("chains", "")) if x.strip()]
    for chain_id in chains:
      out[chain_id] = ChainInfo(chain_id=chain_id, molecule_name=molecule_name, kind=kind)
  return out


def _parse_mmcif_chain_info(path: Path) -> dict[str, ChainInfo]:
  out: dict[str, ChainInfo] = {}
  if path.suffix.lower() not in {".cif", ".mmcif"}:
    return out
  try:
    payload = MMCIF2Dict(str(path))
  except Exception:
    return out

  entity_ids = _as_list(payload.get("_entity.id"))
  entity_descs = _as_list(payload.get("_entity.pdbx_description"))
  entity_name_by_id = {entity_id: _normalize_name(desc) for entity_id, desc in zip(entity_ids, entity_descs, strict=False)}

  poly_entity_ids = _as_list(payload.get("_entity_poly.entity_id"))
  strand_ids = _as_list(payload.get("_entity_poly.pdbx_strand_id"))
  for entity_id, strand_field in zip(poly_entity_ids, strand_ids, strict=False):
    molecule_name = entity_name_by_id.get(entity_id, "")
    kind = _classify_chain_kind(molecule_name)
    chains = [x.strip() for x in strand_field.split(",") if x.strip()]
    for chain_id in chains:
      out[chain_id] = ChainInfo(chain_id=chain_id, molecule_name=molecule_name, kind=kind)
  return out


def _load_chain_info(path: Path) -> dict[str, ChainInfo]:
  if path.suffix.lower() in {".cif", ".mmcif"}:
    return _parse_mmcif_chain_info(path)
  return _parse_pdb_chain_info(path)


def _is_protein_residue(residue) -> bool:
  het, _resseq, _icode = residue.get_id()
  return het == " "


def _split_chains(chain_field: str) -> list[str]:
  raw = (chain_field or "").strip()
  if not raw:
    return []
  for sep in (";", ",", " "):
    if sep in raw:
      return [x.strip() for x in raw.split(sep) if x.strip()]
  return [raw]


def _residue_token(chain_id: str, residue) -> str:
  _het, resseq, icode = residue.get_id()
  icode_s = str(icode).strip()
  if icode_s:
    return f"{chain_id}{int(resseq)}{icode_s}"
  return f"{chain_id}{int(resseq)}"


def _collect_chain_atoms(model, chain_id: str):
  chain = model[chain_id]
  atoms = []
  for residue in chain.get_residues():
    if not _is_protein_residue(residue):
      continue
    atoms.extend(list(residue.get_atoms()))
  return atoms


def _count_contacts_for_chain(model, antigen_chains: list[str], candidate_chain: str, cutoff: float) -> int:
  ab_atoms = _collect_chain_atoms(model, candidate_chain)
  if not ab_atoms:
    return 0
  ns = NeighborSearch(ab_atoms)
  cnt = 0
  for ag_chain in antigen_chains:
    for residue in model[ag_chain].get_residues():
      if not _is_protein_residue(residue):
        continue
      touched = False
      for atom in residue.get_atoms():
        neigh = ns.search(atom.get_coord(), cutoff, level="A")
        if neigh:
          touched = True
          break
      if touched:
        cnt += 1
  return cnt


def _list_available_chains(model) -> list[str]:
  return [c.id for c in model.get_chains()]


def _chain_sequence(chain) -> str:
  seq_chars: list[str] = []
  for residue in chain.get_residues():
    if not _is_protein_residue(residue):
      continue
    seq_chars.append(seq1(residue.get_resname(), custom_map={"MSE": "M"}, undef_code="X"))
  return "".join(seq_chars)


def _load_target_sequence(target_path: Path, requested_chain: str) -> str:
  if not target_path.exists():
    return ""
  try:
    structure = _parse_structure(target_path)
    model = next(iter(structure))
  except Exception:
    return ""

  available = _list_available_chains(model)
  if requested_chain and requested_chain in available:
    return _chain_sequence(model[requested_chain])
  if len(available) == 1:
    return _chain_sequence(model[available[0]])
  return ""


def _has_valid_target_structure(target_path: Path, requested_chain: str) -> bool:
  if not target_path.exists():
    return False
  try:
    structure = _parse_structure(target_path)
    model = next(iter(structure))
  except Exception:
    return False
  available = _list_available_chains(model)
  chains = [requested_chain] if requested_chain and requested_chain in available else available
  for chain_id in chains:
    if any(_is_protein_residue(residue) for residue in model[chain_id].get_residues()):
      return True
  return False


def _infer_requested_binder_kind(row: dict[str, str], requested_chain: str, chain_info: dict[str, ChainInfo]) -> str | None:
  try:
    start = int((row.get("cdr_h3_start") or "").strip())
    end = int((row.get("cdr_h3_end") or "").strip())
  except Exception:
    start = None
    end = None

  if start is not None and end is not None:
    if start <= 92 and end <= 100:
      return "light"
    if start >= 97 or end >= 103:
      return "heavy"

  info = chain_info.get(requested_chain)
  if info and info.kind in BINDER_KINDS:
    return info.kind
  return None


def _candidate_binder_chains(available: list[str], chain_info: dict[str, ChainInfo], antigen_chains: list[str]) -> list[str]:
  binder_like = [chain_id for chain_id in available if chain_id not in antigen_chains and chain_info.get(chain_id, ChainInfo(chain_id, "", "other")).kind in BINDER_KINDS]
  if binder_like:
    return binder_like
  return [chain_id for chain_id in available if chain_id not in antigen_chains]


def _candidate_target_chains(available: list[str], chain_info: dict[str, ChainInfo]) -> list[str]:
  non_binder = [chain_id for chain_id in available if chain_info.get(chain_id, ChainInfo(chain_id, "", "other")).kind == "other"]
  if non_binder:
    return non_binder
  return available[:]


def _match_target_candidates_by_sequence(model, candidates: list[str], target_sequence: str) -> list[str]:
  if not target_sequence:
    return candidates
  exact = [chain_id for chain_id in candidates if _chain_sequence(model[chain_id]) == target_sequence]
  if exact:
    return exact
  return candidates


def _restrict_binder_candidates_for_antigens(model, antigen_chains: list[str], candidates: list[str], requested_kind: str | None, chain_info: dict[str, ChainInfo], cutoff: float) -> list[str]:
  if not requested_kind:
    return candidates
  if requested_kind in {"heavy", "light"}:
    typed = [chain_id for chain_id in candidates if chain_info.get(chain_id, ChainInfo(chain_id, "", "other")).kind == requested_kind]
  else:
    typed = [chain_id for chain_id in candidates if chain_info.get(chain_id, ChainInfo(chain_id, "", "other")).kind in BINDER_KINDS]
  if typed and any(_count_contacts_for_chain(model, antigen_chains, chain_id, cutoff=cutoff) > 0 for chain_id in typed):
    return typed
  return candidates


def _score_antibody_candidate(chain_id: str, contact_count: int, requested_chain: str, requested_kind: str | None, chain_info: dict[str, ChainInfo]) -> int:
  score = contact_count * 100
  info = chain_info.get(chain_id)
  requested_info = chain_info.get(requested_chain)
  if requested_chain:
    if chain_id == requested_chain:
      score += 250
    elif requested_info and info and requested_info.molecule_name and requested_info.molecule_name == info.molecule_name:
      score += 140
  if requested_kind:
    if requested_kind in {"heavy", "light"} and info and info.kind == requested_kind:
      score += 100
    elif requested_kind == "binder" and info and info.kind in BINDER_KINDS:
      score += 60
  return score


def _resolve_antibody_for_fixed_antigen(
  model,
  antigen_chains: list[str],
  requested_chain: str,
  requested_kind: str | None,
  chain_info: dict[str, ChainInfo],
  cutoff: float,
) -> tuple[str | None, str]:
  available = _list_available_chains(model)
  candidates = _candidate_binder_chains(available, chain_info, antigen_chains)
  if not candidates:
    return None, "no_candidate_antibody_chain"
  candidates = _restrict_binder_candidates_for_antigens(model, antigen_chains, candidates, requested_kind, chain_info, cutoff)

  best_chain: str | None = None
  best_score: int | None = None
  best_contact = 0
  for chain_id in candidates:
    contact_count = _count_contacts_for_chain(model, antigen_chains, chain_id, cutoff=cutoff)
    if contact_count <= 0:
      continue
    score = _score_antibody_candidate(chain_id, contact_count, requested_chain, requested_kind, chain_info)
    if best_score is None or score > best_score:
      best_chain = chain_id
      best_score = score
      best_contact = contact_count

  if best_chain is None:
    return None, f"antibody_chain_not_found:need={requested_chain};available={','.join(available)}"
  if best_chain == requested_chain:
    return best_chain, ""
  return best_chain, f"fallback_antibody_chain:{requested_chain}->{best_chain};contact_residues={best_contact}"


def _score_chain_pair(
  antigen_chain: str,
  antibody_chain: str,
  contact_count: int,
  requested_antigen_chains: list[str],
  requested_antibody_chain: str,
  requested_binder_kind: str | None,
  matched_target_chains: list[str],
  chain_info: dict[str, ChainInfo],
) -> int:
  score = _score_antibody_candidate(antibody_chain, contact_count, requested_antibody_chain, requested_binder_kind, chain_info)
  if antigen_chain in requested_antigen_chains:
    score += 200
  if matched_target_chains and antigen_chain in matched_target_chains:
    score += 150
  return score


def _resolve_chain_pair(
  model,
  requested_antigen_chains: list[str],
  requested_antibody_chain: str,
  requested_binder_kind: str | None,
  target_sequence: str,
  chain_info: dict[str, ChainInfo],
  cutoff: float,
) -> tuple[list[str], str | None, str]:
  available = _list_available_chains(model)
  target_candidates = _candidate_target_chains(available, chain_info)
  if not target_candidates:
    return [], None, "no_candidate_antigen_chain"
  target_candidates = _match_target_candidates_by_sequence(model, target_candidates, target_sequence)
  binder_candidates = _candidate_binder_chains(available, chain_info, [])
  if not binder_candidates:
    return [], None, "no_candidate_antibody_chain"
  binder_candidates = _restrict_binder_candidates_for_antigens(model, target_candidates, binder_candidates, requested_binder_kind, chain_info, cutoff)

  best_antigen: str | None = None
  best_antibody: str | None = None
  best_score: int | None = None
  best_contact = 0
  for antigen_chain in target_candidates:
    for antibody_chain in binder_candidates:
      if antigen_chain == antibody_chain:
        continue
      contact_count = _count_contacts_for_chain(model, [antigen_chain], antibody_chain, cutoff=cutoff)
      if contact_count <= 0:
        continue
      score = _score_chain_pair(
        antigen_chain=antigen_chain,
        antibody_chain=antibody_chain,
        contact_count=contact_count,
        requested_antigen_chains=requested_antigen_chains,
        requested_antibody_chain=requested_antibody_chain,
        requested_binder_kind=requested_binder_kind,
        matched_target_chains=target_candidates if target_sequence else [],
        chain_info=chain_info,
      )
      if best_score is None or score > best_score:
        best_antigen = antigen_chain
        best_antibody = antibody_chain
        best_score = score
        best_contact = contact_count

  if best_antigen is None or best_antibody is None:
    return [], None, f"antigen_chain_not_found:need={','.join(requested_antigen_chains)};available={','.join(available)}"

  notes: list[str] = []
  requested_antigen = requested_antigen_chains[0] if requested_antigen_chains else ""
  if requested_antigen and best_antigen != requested_antigen:
    notes.append(f"fallback_antigen_chain:{requested_antigen}->{best_antigen}")
  if requested_antibody_chain and best_antibody != requested_antibody_chain:
    notes.append(f"fallback_antibody_chain:{requested_antibody_chain}->{best_antibody}")
  notes.append(f"contact_residues={best_contact}")
  return [best_antigen], best_antibody, ";".join(notes)


def _write_target_chain_from_complex(ref_path: Path, source_chain_id: str, out_path: Path, out_chain_id: str) -> None:
  structure = _parse_structure(ref_path)
  model = next(iter(structure))
  chain = deepcopy(model[source_chain_id])
  chain.id = out_chain_id

  target_structure = Structure("target")
  target_model = Model(0)
  target_model.add(chain)
  target_structure.add(target_model)

  out_path.parent.mkdir(parents=True, exist_ok=True)
  io = PDBIO()
  io.set_structure(target_structure)
  io.save(str(out_path))


def _extract_hotspots(model, antigen_chains: list[str], antibody_chain: str, output_chain_map: dict[str, str], cutoff: float) -> list[str]:
  ab_atoms = _collect_chain_atoms(model, antibody_chain)
  if not ab_atoms:
    return []
  ns = NeighborSearch(ab_atoms)
  hotspots: list[str] = []

  for ag_chain in antigen_chains:
    for residue in model[ag_chain].get_residues():
      if not _is_protein_residue(residue):
        continue
      hit = False
      for atom in residue.get_atoms():
        neigh = ns.search(atom.get_coord(), cutoff, level="A")
        if neigh:
          hit = True
          break
      if hit:
        output_chain = output_chain_map.get(ag_chain, ag_chain)
        hotspots.append(_residue_token(output_chain, residue))
  # sorted unique while keeping numeric-ish order by insertion in structure
  seen = set()
  out = []
  for x in hotspots:
    if x not in seen:
      seen.add(x)
      out.append(x)
  return out


def main() -> int:
  parser = argparse.ArgumentParser(description="Extract antigen interface hotspots from reference complexes.")
  parser.add_argument(
    "--dataset-csv",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "prepared" / "dataset_index_ready.csv",
    help="Dataset CSV with sample_id/reference_complex_path/antigen_chain/antibody_chain.",
  )
  parser.add_argument(
    "--out-dir",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "prepared" / "epitopes",
    help="Output directory for per-sample epitope JSON files.",
  )
  parser.add_argument(
    "--cutoff",
    type=float,
    default=5.0,
    help="Heavy-atom distance cutoff in Angstrom for interface hotspot extraction.",
  )
  parser.add_argument(
    "--repair-target-pdbs",
    action="store_true",
    help="When reference_structure_path is missing/empty, regenerate the target chain from reference_complex_path.",
  )
  args = parser.parse_args()

  with args.dataset_csv.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))

  args.out_dir.mkdir(parents=True, exist_ok=True)
  summary_rows: list[dict[str, str | int | float]] = []

  ok = 0
  fail = 0
  repaired_targets = 0

  for row in rows:
    sample_id = (row.get("sample_id") or "").strip()
    if not sample_id:
      continue

    ref_path = Path((row.get("reference_complex_path") or row.get("reference_structure_path") or "").strip())
    target_path = Path((row.get("reference_structure_path") or "").strip()) if (row.get("reference_structure_path") or "").strip() else None
    antigen_chains = _split_chains((row.get("antigen_chain") or "").strip())
    ab_chain_req = (row.get("antibody_chain") or "").strip()

    result = SampleResult(
      sample_id=sample_id,
      status="fail",
      note="",
      antigen_chain_used="",
      antibody_chain_used="",
      hotspot_count=0,
      hotspots=[],
      target_pdb_repaired=False,
    )

    if not ref_path.exists():
      result.note = "reference_complex_not_found"
    elif not antigen_chains:
      result.note = "missing_antigen_chain"
    else:
      try:
        structure = _parse_structure(ref_path)
        model = next(iter(structure))
        available = _list_available_chains(model)
        chain_info = _load_chain_info(ref_path)
        requested_binder_kind = _infer_requested_binder_kind(row, ab_chain_req, chain_info)
        target_sequence = _load_target_sequence(target_path, antigen_chains[0]) if target_path is not None and len(antigen_chains) == 1 else ""

        if all(chain_id in available for chain_id in antigen_chains):
          antigen_chains_used = antigen_chains
          ab_chain_used, chain_note = _resolve_antibody_for_fixed_antigen(
            model,
            antigen_chains=antigen_chains_used,
            requested_chain=ab_chain_req,
            requested_kind=requested_binder_kind,
            chain_info=chain_info,
            cutoff=args.cutoff,
          )
        else:
          antigen_chains_used, ab_chain_used, chain_note = _resolve_chain_pair(
            model,
            requested_antigen_chains=antigen_chains,
            requested_antibody_chain=ab_chain_req,
            requested_binder_kind=requested_binder_kind,
            target_sequence=target_sequence,
            chain_info=chain_info,
            cutoff=args.cutoff,
          )

        if not antigen_chains_used:
          result.note = chain_note or f"antigen_chain_not_found:need={','.join(antigen_chains)};available={','.join(available)}"
        elif not ab_chain_used:
          result.note = chain_note or "antibody_chain_not_resolved"
        else:
          output_chain_map = {resolved: antigen_chains[idx] for idx, resolved in enumerate(antigen_chains_used[:len(antigen_chains)])}
          if len(antigen_chains) == 1 and antigen_chains_used:
            output_chain_map = {antigen_chains_used[0]: antigen_chains[0]}

          should_repair_target = (
            args.repair_target_pdbs
            and target_path is not None
            and len(antigen_chains) == 1
            and antigen_chains_used
            and (
              not _has_valid_target_structure(target_path, antigen_chains[0])
              or antigen_chains[0] not in available
            )
          )
          if should_repair_target:
            _write_target_chain_from_complex(ref_path, antigen_chains_used[0], target_path, antigen_chains[0])
            result.target_pdb_repaired = True

          if not ab_chain_used:
            result.note = chain_note or "antibody_chain_not_resolved"
          else:
            hotspots = _extract_hotspots(model, antigen_chains_used, ab_chain_used, output_chain_map, cutoff=args.cutoff)
            result.status = "ok"
            result.note = chain_note
            result.antigen_chain_used = ",".join(antigen_chains_used)
            result.antibody_chain_used = ab_chain_used
            result.hotspots = hotspots
            result.hotspot_count = len(hotspots)
      except StopIteration:
        result.note = "empty_or_invalid_structure"
      except Exception as e:
        result.note = f"parse_failed:{type(e).__name__}"

    out_json = args.out_dir / f"{sample_id}.json"
    payload = {
      "sample_id": sample_id,
      "pdb_id": (row.get("pdb_id") or "").strip(),
      "reference_complex_path": str(ref_path),
      "reference_structure_path": str(target_path) if target_path is not None else "",
      "antigen_chain": (row.get("antigen_chain") or "").strip(),
      "antigen_chain_requested": (row.get("antigen_chain") or "").strip(),
      "antigen_chain_used": result.antigen_chain_used,
      "antibody_chain_requested": ab_chain_req,
      "antibody_chain_used": result.antibody_chain_used,
      "cutoff_angstrom": args.cutoff,
      "status": result.status,
      "note": result.note,
      "hotspot_count": result.hotspot_count,
      "hotspot_residues": result.hotspots,
      "hotspot_string": ",".join(result.hotspots),
      "target_pdb_repaired": result.target_pdb_repaired,
    }
    out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    summary_rows.append(
      {
        "sample_id": sample_id,
        "pdb_id": (row.get("pdb_id") or "").strip(),
        "status": result.status,
        "note": result.note,
        "antigen_chain_requested": (row.get("antigen_chain") or "").strip(),
        "antigen_chain_used": result.antigen_chain_used,
        "antibody_chain_requested": ab_chain_req,
        "antibody_chain_used": result.antibody_chain_used,
        "hotspot_count": result.hotspot_count,
        "hotspot_string": ",".join(result.hotspots),
        "target_pdb_repaired": str(result.target_pdb_repaired),
      }
    )

    if result.status == "ok":
      ok += 1
    else:
      fail += 1
    if result.target_pdb_repaired:
      repaired_targets += 1

  summary_csv = args.out_dir / "epitopes_summary.csv"
  with summary_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(
      f,
      fieldnames=[
        "sample_id",
        "pdb_id",
        "status",
        "note",
        "antigen_chain_requested",
        "antigen_chain_used",
        "antibody_chain_requested",
        "antibody_chain_used",
        "hotspot_count",
        "hotspot_string",
        "target_pdb_repaired",
      ],
    )
    writer.writeheader()
    writer.writerows(summary_rows)

  print(f"[OK] epitope files written: {args.out_dir}")
  print(f"[INFO] summary: {summary_csv}")
  print(f"[INFO] total={len(summary_rows)} ok={ok} fail={fail}")
  print(f"[INFO] repaired_target_pdbs={repaired_targets}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

