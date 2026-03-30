#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import sys

from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1


# IMGT heavy-chain CDR3 interval (inclusive)
H3_RANGES = {
  "imgt": (105, 117),
}


@dataclass
class ChainSeq:
  sequence: str
  residue_map: list[tuple[int, str]]  # sequence index -> (resseq, icode)


def _parse_structure(struct_path: Path):
  if struct_path.suffix.lower() in {".cif", ".mmcif"}:
    parser = MMCIFParser(QUIET=True)
  else:
    parser = PDBParser(QUIET=True)
  return parser.get_structure("ref", str(struct_path))


def _load_chain_sequence(struct_path: Path, chain_id: str) -> ChainSeq:
  structure = _parse_structure(struct_path)
  model = next(iter(structure))
  chain = model[chain_id]

  seq_chars: list[str] = []
  residue_map: list[tuple[int, str]] = []
  for residue in chain.get_residues():
    het, resseq, icode = residue.get_id()
    if het != " ":
      continue
    aa = seq1(residue.get_resname(), custom_map={"MSE": "M"}, undef_code="X")
    seq_chars.append(aa)
    residue_map.append((int(resseq), str(icode) if icode else " "))
  return ChainSeq(sequence="".join(seq_chars), residue_map=residue_map)


def _list_chain_ids(struct_path: Path) -> list[str]:
  structure = _parse_structure(struct_path)
  model = next(iter(structure))
  return [c.id for c in model.get_chains()]


def _resolve_antibody_chain(requested_chain: str, available_chains: list[str]) -> tuple[str | None, str]:
  if requested_chain in available_chains:
    return requested_chain, ""
  # Common PDB convention in antibody complexes.
  for cand in ("H", "V", "L"):
    if cand in available_chains:
      return cand, f"fallback_chain:{requested_chain}->{cand}"
  return None, ""


def _find_position_aa_pairs(obj: Any) -> list[tuple[tuple[int, str], str]] | None:
  """
  递归从 ANARCI 返回结构中提取形如 [((105, ' '), 'A'), ...] 的编号列表。
  """
  if isinstance(obj, list):
    ok = True
    pairs: list[tuple[tuple[int, str], str]] = []
    for x in obj:
      if not (isinstance(x, tuple) and len(x) == 2):
        ok = False
        break
      pos, aa = x
      if not (
        isinstance(pos, tuple)
        and len(pos) >= 2
        and isinstance(pos[0], int)
        and isinstance(pos[1], str)
        and isinstance(aa, str)
      ):
        ok = False
        break
      pairs.append(((int(pos[0]), str(pos[1])), aa))
    if ok:
      return pairs

    for x in obj:
      got = _find_position_aa_pairs(x)
      if got is not None:
        return got

  elif isinstance(obj, tuple):
    for x in obj:
      got = _find_position_aa_pairs(x)
      if got is not None:
        return got

  return None


def _anarci_numbering_pairs(seq: str, scheme: str) -> list[tuple[tuple[int, str], str]] | None:
  try:
    from anarci import anarci as anarci_fn
  except Exception:
    # Fallback: local vendored dependency inside repository.
    repo_root = Path(__file__).resolve().parents[1]
    local_dep = repo_root / ".pydeps_anarci"
    if local_dep.exists():
      sys.path.insert(0, str(local_dep))
      try:
        from anarci import anarci as anarci_fn
      except Exception as e:
        raise RuntimeError(
          "未检测到可用 ANARCI。可尝试 `python -m pip install --target ./.pydeps_anarci --no-deps anarci`。"
        ) from e
    else:
      raise RuntimeError(
        "未检测到 ANARCI Python 包。请先安装（例如 `pip install anarci` 或按 ANARCI 官方方式安装）。"
      )

  # Typical return: (numbering, alignment_details, hit_tables)
  result = anarci_fn([("query", seq)], scheme=scheme)
  if not result or len(result) < 1:
    return None
  numbering = result[0]
  return _find_position_aa_pairs(numbering)


def _infer_h3_start_end(
  seq: str,
  residue_map: list[tuple[int, str]],
  *,
  scheme: str,
  skip_insertions: bool,
) -> tuple[int | None, int | None, str]:
  if scheme not in H3_RANGES:
    return None, None, f"unsupported_scheme:{scheme}"

  pairs = _anarci_numbering_pairs(seq, scheme=scheme)
  if not pairs:
    return None, None, "anarci_no_numbering"

  h3_lo, h3_hi = H3_RANGES[scheme]
  seq_idx = -1
  matched: list[tuple[int, str]] = []

  for (num, ins), aa in pairs:
    if aa != "-":
      seq_idx += 1
    if aa == "-":
      continue
    if not (h3_lo <= num <= h3_hi):
      continue
    if skip_insertions and ins.strip():
      return None, None, "h3_has_anarci_insertion"
    if seq_idx < 0 or seq_idx >= len(residue_map):
      continue
    matched.append(residue_map[seq_idx])

  if not matched:
    return None, None, "h3_not_found"

  if skip_insertions and any(icode.strip() for _resseq, icode in matched):
    return None, None, "h3_has_pdb_insertion"

  start = min(x[0] for x in matched)
  end = max(x[0] for x in matched)
  return start, end, "ok"


def main() -> int:
  parser = argparse.ArgumentParser(description="使用 ANARCI 自动回填 dataset_index.csv 的 CDR H3 起止")
  parser.add_argument("--input-csv", type=Path, required=True, help="输入 dataset_index.csv")
  parser.add_argument("--output-csv", type=Path, required=True, help="输出 CSV（会保留原列并更新 cdr_h3_*）")
  parser.add_argument("--scheme", default="imgt", choices=["imgt"], help="抗体编号方案")
  parser.add_argument(
    "--allow-insertions",
    action="store_true",
    help="允许带插入号的样本写入 H3（默认不允许）",
  )
  parser.add_argument(
    "--force",
    action="store_true",
    help="即使 cdr_h3_start/end 已有值，也强制覆盖",
  )
  args = parser.parse_args()

  with args.input_csv.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))
    fieldnames = list(rows[0].keys()) if rows else []

  for k in ("cdr_h3_start", "cdr_h3_end"):
    if k not in fieldnames:
      fieldnames.append(k)

  # 仅在输出里添加状态列，方便审计；不会影响原有流程。
  for k in ("cdr_h3_status", "cdr_h3_note"):
    if k not in fieldnames:
      fieldnames.append(k)

  total = len(rows)
  updated = 0
  skipped = 0
  failed = 0

  for row in rows:
    if not args.force and (row.get("cdr_h3_start") or row.get("cdr_h3_end")):
      row["cdr_h3_status"] = "skip_existing"
      row["cdr_h3_note"] = ""
      skipped += 1
      continue

    # 推荐：reference_complex_path 放完整复合物；reference_structure_path 常用于模型输入（可能仅抗原）
    ref_path_raw = (row.get("reference_complex_path") or "").strip()
    if not ref_path_raw:
      ref_path_raw = (row.get("reference_structure_path") or "").strip()
    chain_id = (row.get("antibody_chain") or "").strip()
    if not ref_path_raw or not chain_id:
      row["cdr_h3_status"] = "fail"
      row["cdr_h3_note"] = "missing_reference_or_chain"
      failed += 1
      continue

    ref_path = Path(ref_path_raw)
    if not ref_path.exists():
      row["cdr_h3_status"] = "fail"
      row["cdr_h3_note"] = "reference_not_found"
      failed += 1
      continue

    try:
      chain_ids = _list_chain_ids(ref_path)
    except StopIteration:
      row["cdr_h3_status"] = "fail"
      row["cdr_h3_note"] = "empty_or_invalid_structure"
      failed += 1
      continue
    except Exception as e:
      row["cdr_h3_status"] = "fail"
      row["cdr_h3_note"] = f"parse_chain_failed:{type(e).__name__}"
      failed += 1
      continue

    resolved_chain, chain_note = _resolve_antibody_chain(chain_id, chain_ids)
    if resolved_chain is None:
      row["cdr_h3_status"] = "fail"
      row["cdr_h3_note"] = f"antibody_chain_not_found:need={chain_id};available={','.join(chain_ids)}"
      failed += 1
      continue

    try:
      chain_seq = _load_chain_sequence(ref_path, resolved_chain)
    except KeyError:
      row["cdr_h3_status"] = "fail"
      row["cdr_h3_note"] = f"antibody_chain_not_found:need={chain_id};available={','.join(chain_ids)}"
      failed += 1
      continue
    except StopIteration:
      row["cdr_h3_status"] = "fail"
      row["cdr_h3_note"] = "empty_or_invalid_structure"
      failed += 1
      continue
    except Exception as e:
      row["cdr_h3_status"] = "fail"
      row["cdr_h3_note"] = f"parse_chain_failed:{type(e).__name__}"
      failed += 1
      continue

    try:
      h3_start, h3_end, note = _infer_h3_start_end(
        chain_seq.sequence,
        chain_seq.residue_map,
        scheme=args.scheme,
        skip_insertions=(not args.allow_insertions),
      )
    except RuntimeError as e:
      raise SystemExit(str(e))
    except FileNotFoundError:
      row["cdr_h3_status"] = "fail"
      row["cdr_h3_note"] = "anarci_hmmscan_missing"
      failed += 1
      continue
    except Exception as e:
      row["cdr_h3_status"] = "fail"
      row["cdr_h3_note"] = f"anarci_error:{type(e).__name__}"
      failed += 1
      continue

    if h3_start is None or h3_end is None:
      row["cdr_h3_status"] = "skip"
      row["cdr_h3_note"] = f"{chain_note};{note}".strip(";")
      skipped += 1
      continue

    row["cdr_h3_start"] = str(h3_start)
    row["cdr_h3_end"] = str(h3_end)
    row["cdr_h3_status"] = "ok"
    row["cdr_h3_note"] = chain_note
    updated += 1

  args.output_csv.parent.mkdir(parents=True, exist_ok=True)
  with args.output_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

  print(
    f"[OK] 写出: {args.output_csv}\n"
    f"[INFO] total={total} updated={updated} skipped={skipped} failed={failed}"
  )
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

