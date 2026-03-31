#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser


@dataclass
class SampleResult:
  sample_id: str
  status: str
  note: str
  antibody_chain_used: str
  hotspot_count: int
  hotspots: list[str]


def _parse_structure(path: Path):
  if path.suffix.lower() in {".cif", ".mmcif"}:
    parser = MMCIFParser(QUIET=True)
  else:
    parser = PDBParser(QUIET=True)
  return parser.get_structure("complex", str(path))


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


def _resolve_antibody_chain(model, antigen_chains: list[str], requested_chain: str, cutoff: float) -> tuple[str | None, str]:
  available = [c.id for c in model.get_chains()]
  if requested_chain and requested_chain in available and requested_chain not in antigen_chains:
    return requested_chain, ""

  # Fallback: choose the non-antigen chain with highest contact count.
  candidates = [c for c in available if c not in antigen_chains]
  if not candidates:
    return None, "no_candidate_antibody_chain"

  scored = []
  for cand in candidates:
    score = _count_contacts_for_chain(model, antigen_chains, cand, cutoff=cutoff)
    scored.append((score, cand))
  scored.sort(reverse=True)
  best_score, best_chain = scored[0]
  if best_score <= 0:
    return None, f"antibody_chain_not_found:need={requested_chain};available={','.join(available)}"
  return best_chain, f"fallback_antibody_chain:{requested_chain}->{best_chain}"


def _extract_hotspots(model, antigen_chains: list[str], antibody_chain: str, cutoff: float) -> list[str]:
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
        hotspots.append(_residue_token(ag_chain, residue))
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
    default=Path(__file__).resolve().parents[2] / "inputs" / "antibody_datasets" / "dataset_index_ready.csv",
    help="Dataset CSV with sample_id/reference_complex_path/antigen_chain/antibody_chain.",
  )
  parser.add_argument(
    "--out-dir",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "raw" / "epitopes",
    help="Output directory for per-sample epitope JSON files.",
  )
  parser.add_argument(
    "--cutoff",
    type=float,
    default=5.0,
    help="Heavy-atom distance cutoff in Angstrom for interface hotspot extraction.",
  )
  args = parser.parse_args()

  with args.dataset_csv.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))

  args.out_dir.mkdir(parents=True, exist_ok=True)
  summary_rows: list[dict[str, str | int | float]] = []

  ok = 0
  fail = 0

  for row in rows:
    sample_id = (row.get("sample_id") or "").strip()
    if not sample_id:
      continue

    ref_path = Path((row.get("reference_complex_path") or row.get("reference_structure_path") or "").strip())
    antigen_chains = _split_chains((row.get("antigen_chain") or "").strip())
    ab_chain_req = (row.get("antibody_chain") or "").strip()

    result = SampleResult(
      sample_id=sample_id,
      status="fail",
      note="",
      antibody_chain_used="",
      hotspot_count=0,
      hotspots=[],
    )

    if not ref_path.exists():
      result.note = "reference_complex_not_found"
    elif not antigen_chains:
      result.note = "missing_antigen_chain"
    else:
      try:
        structure = _parse_structure(ref_path)
        model = next(iter(structure))
        available = [c.id for c in model.get_chains()]
        if any(c not in available for c in antigen_chains):
          result.note = f"antigen_chain_not_found:need={','.join(antigen_chains)};available={','.join(available)}"
        else:
          ab_chain_used, chain_note = _resolve_antibody_chain(model, antigen_chains, ab_chain_req, cutoff=args.cutoff)
          if not ab_chain_used:
            result.note = chain_note or "antibody_chain_not_resolved"
          else:
            hotspots = _extract_hotspots(model, antigen_chains, ab_chain_used, cutoff=args.cutoff)
            result.status = "ok"
            result.note = chain_note
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
      "antigen_chain": (row.get("antigen_chain") or "").strip(),
      "antibody_chain_requested": ab_chain_req,
      "antibody_chain_used": result.antibody_chain_used,
      "cutoff_angstrom": args.cutoff,
      "status": result.status,
      "note": result.note,
      "hotspot_count": result.hotspot_count,
      "hotspot_residues": result.hotspots,
      "hotspot_string": ",".join(result.hotspots),
    }
    out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    summary_rows.append(
      {
        "sample_id": sample_id,
        "pdb_id": (row.get("pdb_id") or "").strip(),
        "status": result.status,
        "note": result.note,
        "antibody_chain_requested": ab_chain_req,
        "antibody_chain_used": result.antibody_chain_used,
        "hotspot_count": result.hotspot_count,
        "hotspot_string": ",".join(result.hotspots),
      }
    )

    if result.status == "ok":
      ok += 1
    else:
      fail += 1

  summary_csv = args.out_dir / "epitopes_summary.csv"
  with summary_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(
      f,
      fieldnames=[
        "sample_id",
        "pdb_id",
        "status",
        "note",
        "antibody_chain_requested",
        "antibody_chain_used",
        "hotspot_count",
        "hotspot_string",
      ],
    )
    writer.writeheader()
    writer.writerows(summary_rows)

  print(f"[OK] epitope files written: {args.out_dir}")
  print(f"[INFO] summary: {summary_csv}")
  print(f"[INFO] total={len(summary_rows)} ok={ok} fail={fail}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

