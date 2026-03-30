from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select, Superimposer


_BACKBONE_ATOMS_DEFAULT: tuple[str, ...] = ("N", "CA", "C")


@dataclass(frozen=True)
class RMSDResult:
  rmsd: float
  atom_count: int


def _parse_structure(struct_path: Path):
  struct_path = Path(struct_path)
  if not struct_path.exists():
    raise FileNotFoundError(str(struct_path))

  if struct_path.suffix.lower() in {".cif", ".mmcif"}:
    parser = MMCIFParser(QUIET=True)
  else:
    parser = PDBParser(QUIET=True)
  return parser.get_structure("s", str(struct_path))


def _iter_model_chains(structure):
  # Most predictions are single-model; keep it simple and deterministic.
  model = next(iter(structure))
  return model, list(model.get_chains())


def _chain_id_exists(chains, chain_id: str) -> bool:
  return any(c.id == chain_id for c in chains)


def _pick_best_chain_for_antigen(model_chains: list, ref_keys: set[tuple[int, str]], desired_chain_id: str | None):
  # Preferred: requested chain id.
  if desired_chain_id and _chain_id_exists(model_chains, desired_chain_id):
    for c in model_chains:
      if c.id == desired_chain_id:
        return c

  # Fallback: choose chain with maximum overlap by residue keys.
  best = None
  best_score = -1
  for chain in model_chains:
    keys = _collect_residue_keys(chain)
    score = len(keys & ref_keys)
    if score > best_score:
      best_score = score
      best = chain
  return best


def _collect_residue_keys(chain) -> set[tuple[int, str]]:
  # Residue id in Biopython: (hetfield, resseq, icode)
  keys: set[tuple[int, str]] = set()
  for residue in chain.get_residues():
    het, resseq, icode = residue.get_id()
    if het != " ":
      continue
    keys.add((int(resseq), str(icode) if icode != "" else " "))
  return keys


def _select_chain_for_cdr(model_chains, desired_chain_id: str | None, cdr_start: int, cdr_end: int):
  if desired_chain_id and _chain_id_exists(model_chains, desired_chain_id):
    for c in model_chains:
      if c.id == desired_chain_id:
        return c

  # Fallback: pick chain with most residues in [start, end].
  best = None
  best_score = -1
  for chain in model_chains:
    score = 0
    for residue in chain.get_residues():
      het, resseq, _icode = residue.get_id()
      if het != " ":
        continue
      resseq = int(resseq)
      if cdr_start <= resseq <= cdr_end:
        if "CA" in residue:
          score += 1
    if score > best_score:
      best_score = score
      best = chain
  return best


def _residue_keys_in_range(chain, start: int, end: int) -> list[tuple[int, str]]:
  keys: list[tuple[int, str]] = []
  for residue in chain.get_residues():
    het, resseq, icode = residue.get_id()
    if het != " ":
      continue
    resseq = int(resseq)
    if start <= resseq <= end:
      keys.append((resseq, str(icode) if icode != "" else " "))
  keys = sorted(set(keys), key=lambda x: (x[0], x[1]))
  return keys


def _residue_map_by_key(chain) -> dict[tuple[int, str], object]:
  out: dict[tuple[int, str], object] = {}
  for residue in chain.get_residues():
    het, resseq, icode = residue.get_id()
    if het != " ":
      continue
    out[(int(resseq), str(icode) if icode != "" else " ")] = residue
  return out


def _compute_superimposition_on_antigen(
  pred_antigen_chain,
  ref_antigen_chain,
  atom_name: str = "CA",
):
  pred_map = _residue_map_by_key(pred_antigen_chain)
  ref_keys = _residue_keys_in_range(ref_antigen_chain, start=-10**9, end=10**9)
  ref_map = _residue_map_by_key(ref_antigen_chain)

  fixed_atoms = []
  moving_atoms = []
  for key in sorted(ref_keys, key=lambda x: (x[0], x[1])):
    ref_res = ref_map.get(key)
    pred_res = pred_map.get(key)
    if ref_res is None or pred_res is None:
      continue
    if atom_name not in ref_res or atom_name not in pred_res:
      continue
    fixed_atoms.append(ref_res[atom_name])
    moving_atoms.append(pred_res[atom_name])

  if len(fixed_atoms) < 3:
    return None

  sup = Superimposer()
  sup.set_atoms(fixed_atoms, moving_atoms)
  rot, tran = sup.rotran  # rot: 3x3, tran: (3,)
  return rot, tran


def _extract_cdr_atom_pairs(
  pred_ab_chain,
  ref_ab_chain,
  cdr_start: int,
  cdr_end: int,
  atom_names: Sequence[str] = _BACKBONE_ATOMS_DEFAULT,
):
  pred_map = _residue_map_by_key(pred_ab_chain)
  ref_keys = _residue_keys_in_range(ref_ab_chain, cdr_start, cdr_end)

  pred_coords: list[np.ndarray] = []
  ref_coords: list[np.ndarray] = []
  for key in ref_keys:
    ref_res = _residue_map_by_key(ref_ab_chain).get(key)
    pred_res = pred_map.get(key)
    if ref_res is None or pred_res is None:
      continue
    for atom_name in atom_names:
      if atom_name in ref_res and atom_name in pred_res:
        pred_coords.append(np.asarray(pred_res[atom_name].get_coord(), dtype=float))
        ref_coords.append(np.asarray(ref_res[atom_name].get_coord(), dtype=float))

  if not pred_coords:
    return np.zeros((0, 3), dtype=float), np.zeros((0, 3), dtype=float), 0
  return np.vstack(pred_coords), np.vstack(ref_coords), len(pred_coords)


def compute_antigen_aligned_cdr_rmsd(
  pred_path: Path,
  ref_path: Path,
  *,
  antigen_chain_id: str,
  antibody_chain_id: str,
  cdr_start: int,
  cdr_end: int,
  cdr_atom_names: Sequence[str] = _BACKBONE_ATOMS_DEFAULT,
) -> Optional[RMSDResult]:
  """
  抗原CA对齐（固定参考复合物抗原），然后在同一刚体变换坐标系中计算抗体CDR的RMSD。

  数据索引建议保证：
  - antigen_chain_id 在 ref 与 pred 中大体对应同一条链（若不一致，会尝试按残基编号重叠自动回退）。
  - antibody_chain_id 在 ref 与 pred 中大体对应同一条链（若不一致，会按CDR区间残基数量自动回退）。
  - cdr_start/cdr_end 使用 ref 的残基编号；pred 会用相同编号提取并建立原子对应关系。
  """
  if cdr_start is None or cdr_end is None:
    return None
  if cdr_start > cdr_end:
    return None

  pred_struct = _parse_structure(pred_path)
  ref_struct = _parse_structure(ref_path)

  _, pred_chains = _iter_model_chains(pred_struct)
  _, ref_chains = _iter_model_chains(ref_struct)

  ref_antigen_chain = None
  for c in ref_chains:
    if c.id == antigen_chain_id:
      ref_antigen_chain = c
      break
  if ref_antigen_chain is None:
    ref_antigen_chain = ref_chains[0] if ref_chains else None
  if ref_antigen_chain is None:
    return None

  ref_antigen_keys = _collect_residue_keys(ref_antigen_chain)
  pred_antigen_chain = _pick_best_chain_for_antigen(
    pred_chains,
    ref_keys=ref_antigen_keys,
    desired_chain_id=antigen_chain_id,
  )
  if pred_antigen_chain is None:
    return None

  rot_tran = _compute_superimposition_on_antigen(pred_antigen_chain, ref_antigen_chain, atom_name="CA")
  if rot_tran is None:
    return None
  rot, tran = rot_tran

  ref_ab_chain = _select_chain_for_cdr(ref_chains, desired_chain_id=antibody_chain_id, cdr_start=cdr_start, cdr_end=cdr_end)
  pred_ab_chain = _select_chain_for_cdr(pred_chains, desired_chain_id=antibody_chain_id, cdr_start=cdr_start, cdr_end=cdr_end)
  if ref_ab_chain is None or pred_ab_chain is None:
    return None

  pred_atoms, ref_atoms, atom_count = _extract_cdr_atom_pairs(
    pred_ab_chain, ref_ab_chain, cdr_start=cdr_start, cdr_end=cdr_end, atom_names=cdr_atom_names
  )
  if atom_count == 0:
    return None

  pred_aligned = (rot @ pred_atoms.T).T + tran  # (N,3)
  diff = pred_aligned - ref_atoms
  rmsd = math.sqrt(float(np.mean(np.sum(diff * diff, axis=1))))
  return RMSDResult(rmsd=rmsd, atom_count=atom_count)


class _CDRSelect(Select):
  def __init__(self, chain_id: str, start: int, end: int) -> None:
    super().__init__()
    self.chain_id = chain_id
    self.start = start
    self.end = end

  def accept_chain(self, chain) -> int:
    return 1 if chain.id == self.chain_id else 0

  def accept_residue(self, residue) -> int:
    het, resseq, _icode = residue.get_id()
    if het != " ":
      return 0
    resseq = int(resseq)
    return 1 if self.start <= resseq <= self.end else 0


def extract_cdr_backbone_pdb(
  struct_path: Path,
  out_pdb_path: Path,
  *,
  chain_id: str,
  start: int,
  end: int,
) -> bool:
  """
  仅用于人工检查：导出指定链的残基范围（包括全部原子）。
  """
  if start > end:
    return False

  struct = _parse_structure(struct_path)
  io = PDBIO()
  io.set_structure(struct)
  selector = _CDRSelect(chain_id=chain_id, start=start, end=end)
  out_pdb_path = Path(out_pdb_path)
  out_pdb_path.parent.mkdir(parents=True, exist_ok=True)
  io.save(str(out_pdb_path), selector=selector)
  return True

