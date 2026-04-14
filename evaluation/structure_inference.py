from __future__ import annotations

from difflib import SequenceMatcher
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1


def parse_structure(path: Path):
  if path.suffix.lower() in {".cif", ".mmcif"}:
    parser = MMCIFParser(QUIET=True)
  else:
    parser = PDBParser(QUIET=True)
  return parser.get_structure("s", str(path))


def first_model(structure):
  return next(iter(structure))


def protein_chain_sequences(structure) -> dict[str, str]:
  model = first_model(structure)
  sequences: dict[str, str] = {}
  for chain in model.get_chains():
    seq_chars: list[str] = []
    for residue in chain.get_residues():
      if residue.get_id()[0] != " ":
        continue
      seq_chars.append(seq1(residue.get_resname(), custom_map={"MSE": "M"}, undef_code="X"))
    if seq_chars:
      sequences[chain.id] = "".join(seq_chars)
  return sequences


def load_reference_chain_sequence(path: Path | None, preferred_chain_id: str = "") -> str:
  if path is None or not path.exists():
    return ""
  try:
    structure = parse_structure(path)
  except Exception:
    return ""
  sequences = protein_chain_sequences(structure)
  if preferred_chain_id and preferred_chain_id in sequences:
    return sequences[preferred_chain_id]
  if len(sequences) == 1:
    return next(iter(sequences.values()))
  if sequences:
    return max(sequences.values(), key=len)
  return ""


def _sequence_score(candidate_sequence: str, reference_sequence: str) -> tuple[float, int, int]:
  similarity = SequenceMatcher(None, candidate_sequence, reference_sequence).ratio()
  length_delta = abs(len(candidate_sequence) - len(reference_sequence))
  return similarity, -length_delta, len(candidate_sequence)


def _pick_chain_by_sequence(
  chain_sequences: dict[str, str],
  *,
  reference_sequence: str,
  preferred_chain_id: str = "",
  excluded_chain_ids: set[str] | None = None,
  fallback_chain_ids: tuple[str, ...] = (),
  prefer_longest: bool = False,
) -> str | None:
  excluded = excluded_chain_ids or set()
  candidates = [chain_id for chain_id in chain_sequences if chain_id not in excluded]
  if not candidates:
    return None
  if preferred_chain_id and preferred_chain_id in chain_sequences and preferred_chain_id not in excluded:
    return preferred_chain_id
  if reference_sequence:
    ranked = sorted(
      candidates,
      key=lambda chain_id: _sequence_score(chain_sequences[chain_id], reference_sequence),
      reverse=True,
    )
    if ranked:
      return ranked[0]
  for chain_id in fallback_chain_ids:
    if chain_id in chain_sequences and chain_id not in excluded:
      return chain_id
  if prefer_longest:
    return max(candidates, key=lambda chain_id: len(chain_sequences[chain_id]))
  return min(candidates, key=lambda chain_id: len(chain_sequences[chain_id]))


def resolve_structure_chain_ids(
  chain_sequences: dict[str, str],
  *,
  target_sequence: str = "",
  binder_sequence: str = "",
  target_chain_hint: str = "",
  binder_chain_hint: str = "",
) -> tuple[str | None, str | None]:
  target_chain_id = _pick_chain_by_sequence(
    chain_sequences,
    reference_sequence=target_sequence,
    preferred_chain_id=target_chain_hint,
    fallback_chain_ids=("T", "A"),
  )

  binder_chain_id = _pick_chain_by_sequence(
    chain_sequences,
    reference_sequence=binder_sequence,
    preferred_chain_id=binder_chain_hint,
    excluded_chain_ids={target_chain_id} if target_chain_id else set(),
    fallback_chain_ids=("H", "B", "L", "A"),
    prefer_longest=True,
  )

  if target_chain_id is None and binder_chain_id is not None:
    remaining = [chain_id for chain_id in chain_sequences if chain_id != binder_chain_id]
    if len(remaining) == 1:
      target_chain_id = remaining[0]
  if binder_chain_id is None and target_chain_id is not None:
    remaining = [chain_id for chain_id in chain_sequences if chain_id != target_chain_id]
    if remaining:
      binder_chain_id = max(remaining, key=lambda chain_id: len(chain_sequences[chain_id]))
  return target_chain_id, binder_chain_id


def infer_binder_sequence_from_structure(
  *,
  structure_path: Path,
  target_structure_path: Path | None,
  target_chain_hint: str = "",
  binder_chain_hint: str = "",
  binder_sequence_hint: str = "",
) -> tuple[str, str]:
  if not structure_path.exists():
    return "", ""
  try:
    structure = parse_structure(structure_path)
  except Exception:
    return "", ""
  chain_sequences = protein_chain_sequences(structure)
  if not chain_sequences:
    return "", ""
  target_sequence = load_reference_chain_sequence(target_structure_path, preferred_chain_id=target_chain_hint)
  _target_chain_id, binder_chain_id = resolve_structure_chain_ids(
    chain_sequences,
    target_sequence=target_sequence,
    binder_sequence=binder_sequence_hint,
    target_chain_hint=target_chain_hint,
    binder_chain_hint=binder_chain_hint,
  )
  if binder_chain_id is None:
    return "", ""
  return chain_sequences[binder_chain_id], binder_chain_id