#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select

from evaluation.schema import EXTERNAL_CANDIDATE_FIELDS


AA20 = set("ACDEFGHIKLMNPQRSTVWY")
SELECTION_AUDIT_FIELDS = [
  "candidate_id",
  "run_tag",
  "model",
  "sample_id",
  "manifest_row_idx",
  "source_candidate_rank_in_model",
  "selection_status",
  "selection_reason",
  "selection_detail",
  "pred_sequence",
  "pred_seq_len",
  "source_sequence_path",
  "source_structure_path",
]


class _ChainSelect(Select):
  def __init__(self, chain_id: str) -> None:
    super().__init__()
    self.chain_id = chain_id

  def accept_chain(self, chain) -> int:
    return 1 if chain.id == self.chain_id else 0


def _parse_structure(struct_path: Path):
  if struct_path.suffix.lower() in {".cif", ".mmcif"}:
    parser = MMCIFParser(QUIET=True)
  else:
    parser = PDBParser(QUIET=True)
  return parser.get_structure("s", str(struct_path))


def _safe_read_json(path: Path) -> dict:
  if not path.exists():
    return {}
  try:
    return json.loads(path.read_text(encoding="utf-8"))
  except Exception:
    return {}


def _read_fasta_sequence(path: Path) -> str:
  if not path.exists():
    return ""
  parts: list[str] = []
  with path.open("r", encoding="utf-8") as f:
    for line in f:
      text = line.strip()
      if not text or text.startswith(">"):
        continue
      parts.append(text)
  return "".join(parts).strip().upper()


def _resolve_repo_path(repo_root: Path, raw_path: str) -> Path | None:
  text = (raw_path or "").strip()
  if not text:
    return None
  path = Path(text)
  if not path.is_absolute():
    path = repo_root / path
  return path.resolve()


def _to_int(value: str) -> int | None:
  text = (value or "").strip()
  if not text:
    return None
  try:
    return int(float(text))
  except ValueError:
    return None


def _normalize_sequence(seq: str) -> str:
  return "".join((seq or "").strip().upper().split())


def _load_index_by_sample(index_csv: Path) -> dict[str, dict[str, str]]:
  out: dict[str, dict[str, str]] = {}
  if not index_csv.exists():
    return out
  with index_csv.open("r", encoding="utf-8", newline="") as f:
    for row in csv.DictReader(f):
      sample_id = (row.get("sample_id") or "").strip()
      if sample_id:
        out[sample_id] = row
  return out


def _iter_manifest_csvs(native_roots: list[Path]):
  for root_priority, native_root in enumerate(native_roots):
    run_tag = native_root.name
    if not native_root.exists():
      continue
    for model_dir in sorted(path for path in native_root.iterdir() if path.is_dir()):
      manifest_csv = model_dir / "manifest.csv"
      if manifest_csv.exists():
        yield root_priority, run_tag, model_dir.name, manifest_csv


def _extract_target_chain(reference_complex_path: Path, target_chain_id: str, out_path: Path) -> str:
  if not reference_complex_path.exists() or not target_chain_id:
    return ""
  out_path.parent.mkdir(parents=True, exist_ok=True)
  try:
    structure = _parse_structure(reference_complex_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), _ChainSelect(target_chain_id))
    return str(out_path.resolve())
  except Exception:
    return ""


def _compute_sequence_stats(seq: str) -> tuple[str, str, str, str]:
  if not seq:
    return "", "", "", "0"
  valid_ratio = sum(1 for aa in seq if aa in AA20) / len(seq)
  cys_count = seq.count("C")
  return str(len(seq)), str(round(valid_ratio, 6)), str(cys_count), "1"


def _write_judge_input(job_path: Path, row: dict[str, str]) -> str:
  job_path.parent.mkdir(parents=True, exist_ok=True)
  payload = {
    "candidate_id": row["candidate_id"],
    "sample_id": row["sample_id"],
    "model": row["model"],
    "pred_sequence": row["pred_sequence"],
    "target_structure_path": row["target_structure_path"],
    "target_chain_id": row["target_chain_id"],
    "designed_antibody_chain_id": row["designed_antibody_chain_id"],
    "hotspot_string": row["hotspot_string"],
  }
  job_path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
  return str(job_path.resolve())


def _selection_sort_key(row: dict[str, str]) -> tuple[int, int, int, str]:
  return (
    _to_int(row.get("_root_priority", "")) or 10**9,
    _to_int(row.get("source_candidate_rank_in_model", "")) or 10**9,
    _to_int(row.get("manifest_row_idx", "")) or 10**9,
    row.get("candidate_id", "") or "",
  )


def _build_selection_audit_row(row: dict[str, str]) -> dict[str, str]:
  audit = {field: "" for field in SELECTION_AUDIT_FIELDS}
  for field in [
    "candidate_id",
    "run_tag",
    "model",
    "sample_id",
    "manifest_row_idx",
    "source_candidate_rank_in_model",
    "pred_seq_len",
    "source_sequence_path",
    "source_structure_path",
  ]:
    audit[field] = row.get(field, "") or ""
  audit["pred_sequence"] = _normalize_sequence(row.get("pred_sequence", ""))
  return audit


def _write_selection_audit(path: Path, rows: list[dict[str, str]]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=SELECTION_AUDIT_FIELDS)
    writer.writeheader()
    writer.writerows(rows)


def _select_candidate_rows(
  rows: list[dict[str, str]],
  *,
  max_candidates_per_model_sample: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
  grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
  for row in rows:
    key = ((row.get("model") or "").strip(), (row.get("sample_id") or "").strip())
    grouped.setdefault(key, []).append(row)

  selected_rows: list[dict[str, str]] = []
  audit_rows: list[dict[str, str]] = []
  for key in sorted(grouped):
    kept_by_sequence: dict[str, str] = {}
    selected_count = 0
    for row in sorted(grouped[key], key=_selection_sort_key):
      audit = _build_selection_audit_row(row)
      seq = _normalize_sequence(row.get("pred_sequence", ""))
      if not seq:
        audit["selection_status"] = "excluded"
        audit["selection_reason"] = "missing_pred_sequence"
      elif any(aa not in AA20 for aa in seq):
        audit["selection_status"] = "excluded"
        audit["selection_reason"] = "invalid_pred_sequence"
      elif seq in kept_by_sequence:
        audit["selection_status"] = "excluded"
        audit["selection_reason"] = "duplicate_pred_sequence"
        audit["selection_detail"] = kept_by_sequence[seq]
      elif selected_count >= max_candidates_per_model_sample:
        audit["selection_status"] = "excluded"
        audit["selection_reason"] = "beyond_top_k_unique_sequences"
        audit["selection_detail"] = str(max_candidates_per_model_sample)
      else:
        selected_count += 1
        kept_by_sequence[seq] = row.get("candidate_id", "") or ""
        row["candidate_rank_in_model"] = str(selected_count)
        row["pred_sequence"] = seq
        audit["selection_status"] = "selected"
        audit["selection_reason"] = "top_k_unique_sequence"
        selected_rows.append(row)
      audit_rows.append(audit)
  return selected_rows, audit_rows


def build_external_eval_inputs(
  *,
  native_roots: list[Path],
  index_csv: Path,
  out_csv: Path,
  job_dir: Path,
  judge_name: str = "chai1_primary",
  n_requested_seeds: int = 3,
  max_candidates_per_model_sample: int = 10,
  selection_audit_csv: Path | None = None,
) -> Path:
  repo_root = Path(__file__).resolve().parents[2]
  prepared_dir = repo_root / "data" / "prepared"
  index_by_sample = _load_index_by_sample(index_csv)
  target_dir = job_dir / "targets"
  rows: list[dict[str, str]] = []

  for root_priority, run_tag, model, manifest_csv in _iter_manifest_csvs(native_roots):
    with manifest_csv.open("r", encoding="utf-8", newline="") as f:
      reader = csv.DictReader(f)
      for manifest_row_idx, manifest_row in enumerate(reader, start=1):
        sample_id = (manifest_row.get("sample_id") or "").strip()
        manifest_candidate_id = (manifest_row.get("candidate_id") or "").strip()
        manifest_candidate_rank = (
          manifest_row.get("candidate_rank")
          or manifest_row.get("candidate_rank_in_model")
          or manifest_row.get("candidate_rank_in_sample")
          or ""
        ).strip()
        index_row = index_by_sample.get(sample_id, {})
        sequence_path = _resolve_repo_path(repo_root, manifest_row.get("sequence_path") or "")
        structure_path = _resolve_repo_path(repo_root, manifest_row.get("structure_path") or "")
        meta_path = _resolve_repo_path(repo_root, manifest_row.get("meta_path") or "")
        seq = _read_fasta_sequence(sequence_path) if sequence_path else ""
        pred_seq_len, aa_valid_ratio, cys_count, has_source_sequence = _compute_sequence_stats(seq)

        epitope_path = prepared_dir / "epitopes" / f"{sample_id}.json"
        epitope = _safe_read_json(epitope_path)
        antigen_chain_id = str(epitope.get("antigen_chain") or index_row.get("antigen_chain") or "").strip()
        antibody_chain_id = str(epitope.get("antibody_chain_used") or index_row.get("antibody_chain") or "").strip()
        hotspot_string = str(epitope.get("hotspot_string") or "").strip()
        hotspot_count = str(epitope.get("hotspot_count") or 0)

        reference_complex_path = _resolve_repo_path(repo_root, index_row.get("reference_complex_path") or "")
        target_structure_path = ""
        if reference_complex_path and antigen_chain_id:
          target_structure_path = _extract_target_chain(
            reference_complex_path=reference_complex_path,
            target_chain_id=antigen_chain_id,
            out_path=target_dir / f"{sample_id}_antigen.pdb",
          )

        row = {k: "" for k in EXTERNAL_CANDIDATE_FIELDS}
        row["candidate_id"] = manifest_candidate_id or f"{run_tag}:{model}:{sample_id}:{manifest_row_idx}"
        row["run_id"] = row["candidate_id"]
        row["run_tag"] = run_tag
        row["model"] = model
        row["manifest_row_idx"] = str(manifest_row_idx)
        row["source_candidate_rank_in_model"] = manifest_candidate_rank or str(manifest_row_idx)
        row["candidate_rank_in_model"] = row["source_candidate_rank_in_model"]
        row["sample_id"] = sample_id
        row["source_status"] = (manifest_row.get("status") or "").strip()
        row["source_error_summary"] = (manifest_row.get("error_summary") or "").strip()
        row["source_duration_sec"] = (manifest_row.get("duration_sec") or "").strip()
        row["source_meta_path"] = str(meta_path) if meta_path else ""
        row["source_sequence_path"] = str(sequence_path) if sequence_path else ""
        row["source_structure_path"] = str(structure_path) if structure_path else ""
        row["pred_sequence"] = seq
        row["pred_seq_len"] = pred_seq_len
        row["aa_valid_ratio"] = aa_valid_ratio
        row["cys_count"] = cys_count
        row["has_source_sequence"] = has_source_sequence
        row["has_source_structure"] = "1" if structure_path and structure_path.exists() else "0"
        row["reference_complex_path"] = str(reference_complex_path) if reference_complex_path else ""
        row["target_structure_path"] = target_structure_path
        row["target_chain_id"] = antigen_chain_id
        row["designed_antibody_chain_id"] = antibody_chain_id
        row["epitope_path"] = str(epitope_path) if epitope_path.exists() else ""
        row["hotspot_count"] = hotspot_count
        row["hotspot_string"] = hotspot_string
        row["num_missing_residues"] = ""
        row["judge_name"] = judge_name
        row["n_requested_seeds"] = str(n_requested_seeds)
        row["model_native_metadata_path"] = row["source_meta_path"]
        row["judge_input_path"] = ""
        row["_root_priority"] = str(root_priority)
        rows.append(row)

  selected_rows, audit_rows = _select_candidate_rows(
    rows,
    max_candidates_per_model_sample=max(1, max_candidates_per_model_sample),
  )
  for row in selected_rows:
    row.pop("_root_priority", None)
    candidate_slug = row["candidate_id"].replace(":", "__").replace("/", "_")
    row["judge_input_path"] = _write_judge_input(job_dir / "jobs" / f"{candidate_slug}.json", row)

  out_csv.parent.mkdir(parents=True, exist_ok=True)
  with out_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=EXTERNAL_CANDIDATE_FIELDS)
    writer.writeheader()
    writer.writerows(selected_rows)

  audit_path = selection_audit_csv or out_csv.with_name(f"{out_csv.stem}_selection_audit.csv")
  _write_selection_audit(audit_path, audit_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="构建统一外部评估候选输入表与 judge 任务模板")
  parser.add_argument("--native-root", action="append", type=Path, help="native 输出目录（可重复指定）")
  parser.add_argument("--index-csv", type=Path, default=Path("data/prepared/dataset_index_ready.csv"), help="样本索引 CSV")
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_candidates.csv"),
    help="候选表输出路径",
  )
  parser.add_argument(
    "--job-dir",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/judge_inputs"),
    help="judge 输入模板目录",
  )
  parser.add_argument("--judge-name", type=str, default="chai1_primary", help="主 judge 名称")
  parser.add_argument("--n-requested-seeds", type=int, default=3, help="每条候选请求 seeds 数")
  parser.add_argument(
    "--max-candidates-per-model-sample",
    type=int,
    default=10,
    help="每个 model/sample 保留的 top-K 唯一候选序列数",
  )
  parser.add_argument(
    "--selection-audit-csv",
    type=Path,
    default=None,
    help="候选筛选审计表输出路径，默认与候选表同目录",
  )
  args = parser.parse_args()

  native_roots = args.native_root or [Path("outputs/native_predictions_run2")]
  out_csv = build_external_eval_inputs(
    native_roots=native_roots,
    index_csv=args.index_csv,
    out_csv=args.out_csv,
    job_dir=args.job_dir,
    judge_name=args.judge_name,
    n_requested_seeds=max(1, args.n_requested_seeds),
    max_candidates_per_model_sample=max(1, args.max_candidates_per_model_sample),
    selection_audit_csv=args.selection_audit_csv,
  )
  print(f"[OK] 候选表: {out_csv}")
  print(f"[OK] 候选筛选审计表: {args.selection_audit_csv or args.out_csv.with_name(f'{args.out_csv.stem}_selection_audit.csv')}")
  print(f"[OK] judge 输入目录: {args.job_dir}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
