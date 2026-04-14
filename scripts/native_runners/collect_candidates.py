#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import re
from pathlib import Path


MANIFEST_FIELDS = [
  "candidate_id",
  "sample_id",
  "model",
  "candidate_rank",
  "candidate_name",
  "status",
  "sequence_path",
  "structure_path",
  "meta_path",
  "duration_sec",
  "error_summary",
  "source_stage",
  "source_name",
  "source_csv_path",
  "source_row_idx",
  "source_sequence_field",
  "source_structure_field",
  "native_score_name",
  "native_score_value",
  "notes",
]

FASTA_EXTENSIONS = (".fa", ".faa", ".fasta", ".fas", ".seq")
STRUCTURE_EXTENSIONS = (".pdb", ".cif", ".mmcif")


def _print_fieldnames() -> int:
  print(",".join(MANIFEST_FIELDS))
  return 0


def _slugify(value: str) -> str:
  text = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
  return text.strip("._") or "candidate"


def _normalize_sequence(sequence: str) -> str:
  return "".join((sequence or "").split()).upper()


def _to_project_rel(path: Path | None, project_root: Path) -> str:
  if path is None:
    return ""
  try:
    resolved = path.resolve()
  except FileNotFoundError:
    resolved = path
  try:
    return str(resolved.relative_to(project_root.resolve()))
  except Exception:
    return str(resolved)


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
  if not path.exists():
    return []
  with path.open("r", encoding="utf-8", newline="") as f:
    return [{str(k): str(v) for k, v in row.items()} for row in csv.DictReader(f)]


def _read_json(path: Path) -> dict:
  if not path.exists():
    return {}
  try:
    return json.loads(path.read_text(encoding="utf-8"))
  except Exception:
    return {}


def _write_json(path: Path, payload: dict) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _write_fasta(path: Path, header: str, sequence: str) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_text(f">{header}\n{sequence}\n", encoding="utf-8")


def _resolve_candidate_path(raw: str, base_dir: Path, search_roots: list[Path]) -> Path | None:
  raw = (raw or "").strip()
  if not raw:
    return None
  candidate = Path(raw)
  if candidate.is_absolute() and candidate.exists():
    return candidate
  local = (base_dir / candidate).resolve()
  if local.exists():
    return local
  for root in search_roots:
    alt = (root / candidate).resolve()
    if alt.exists():
      return alt
  return None


def _glob_first(root: Path, patterns: list[str]) -> Path | None:
  for pattern in patterns:
    matches = sorted(root.glob(pattern))
    if matches:
      return matches[0]
  return None


def _glob_all(root: Path, patterns: list[str]) -> list[Path]:
  out: list[Path] = []
  seen: set[str] = set()
  for pattern in patterns:
    for match in sorted(root.glob(pattern)):
      key = str(match.resolve()) if match.exists() else str(match)
      if key in seen:
        continue
      seen.add(key)
      out.append(match)
  return out


def _extract_sequence_from_row(row: dict[str, str], preferred_fields: list[str] | None = None) -> tuple[str, str]:
  preferred = [field.lower() for field in (preferred_fields or [])]
  lowered = {str(k).lower(): str(k) for k in row.keys()}
  for key in preferred:
    actual = lowered.get(key)
    if actual is None:
      continue
    sequence = _normalize_sequence(row.get(actual, ""))
    if sequence:
      return sequence, actual

  generic_priority = [
    "sequence",
    "binder_sequence",
    "designed_sequence",
    "designed_chain_sequence",
    "seq",
    "aa",
  ]
  for key in generic_priority:
    actual = lowered.get(key)
    if actual is None:
      continue
    sequence = _normalize_sequence(row.get(actual, ""))
    if sequence:
      return sequence, actual

  for actual_key, value in row.items():
    key = actual_key.strip().lower()
    if "sequence" not in key:
      continue
    sequence = _normalize_sequence(value)
    if sequence:
      return sequence, actual_key
  return "", ""


def _pick_native_score(row: dict[str, str], preferred_fields: list[str]) -> tuple[str, str]:
  lowered = {str(k).lower(): str(k) for k in row.keys()}
  for key in preferred_fields:
    actual = lowered.get(key.lower())
    if actual is None:
      continue
    value = (row.get(actual) or "").strip()
    if value:
      return actual, value
  return "", ""


def _sort_rows(rows: list[dict[str, str]], rank_field_names: list[str], score_field_names: list[str]) -> list[tuple[int, dict[str, str]]]:
  lowered_rows = []
  for idx, row in enumerate(rows, start=1):
    lowered = {str(k).lower(): str(k) for k in row.keys()}
    rank_value = None
    for field in rank_field_names:
      actual = lowered.get(field.lower())
      if actual is None:
        continue
      raw = (row.get(actual) or "").strip()
      if not raw:
        continue
      try:
        rank_value = float(raw)
        break
      except ValueError:
        continue

    score_value = None
    for field in score_field_names:
      actual = lowered.get(field.lower())
      if actual is None:
        continue
      raw = (row.get(actual) or "").strip()
      if not raw:
        continue
      try:
        score_value = float(raw)
        break
      except ValueError:
        continue

    lowered_rows.append((idx, row, rank_value, score_value))

  def sort_key(item: tuple[int, dict[str, str], float | None, float | None]):
    idx, _row, rank_value, score_value = item
    if rank_value is not None:
      return (0, rank_value, idx)
    if score_value is not None:
      return (1, -score_value, idx)
    return (2, idx, idx)

  lowered_rows.sort(key=sort_key)
  return [(idx, row) for idx, row, _rank_value, _score_value in lowered_rows]


def _build_record(
  *,
  project_root: Path,
  sample_output_dir: Path,
  sample_id: str,
  model: str,
  candidate_rank: int,
  candidate_name: str,
  sequence: str,
  structure_path: Path | None,
  duration_sec: str,
  source_stage: str,
  source_name: str,
  source_csv_path: Path | None,
  source_row_idx: int | None,
  source_sequence_field: str,
  source_structure_field: str,
  native_score_name: str,
  native_score_value: str,
  error_summary: str,
  notes: list[str],
  raw_row: dict[str, str] | None,
) -> dict[str, str]:
  candidate_slug = _slugify(f"{candidate_rank:03d}_{source_stage}_{candidate_name}")
  candidate_id = f"{sample_id}:{model}:{candidate_slug}"
  sequence_path = ""
  normalized_sequence = _normalize_sequence(sequence)
  if normalized_sequence:
    fasta_path = sample_output_dir / "candidate_sequences" / f"{candidate_slug}.fasta"
    _write_fasta(fasta_path, candidate_id, normalized_sequence)
    sequence_path = _to_project_rel(fasta_path, project_root)

  meta_payload = {
    "candidate_id": candidate_id,
    "sample_id": sample_id,
    "model": model,
    "candidate_rank": candidate_rank,
    "candidate_name": candidate_name,
    "source_stage": source_stage,
    "source_name": source_name,
    "source_csv_path": str(source_csv_path) if source_csv_path else "",
    "source_row_idx": source_row_idx,
    "source_sequence_field": source_sequence_field,
    "source_structure_field": source_structure_field,
    "native_score_name": native_score_name,
    "native_score_value": native_score_value,
    "sequence_path": sequence_path,
    "structure_path": _to_project_rel(structure_path, project_root),
    "notes": notes,
    "raw_row": raw_row or {},
  }
  meta_path = sample_output_dir / "candidate_meta" / f"{candidate_slug}.json"
  _write_json(meta_path, meta_payload)

  status = "ok" if sequence_path or structure_path else "failed"
  return {
    "candidate_id": candidate_id,
    "sample_id": sample_id,
    "model": model,
    "candidate_rank": str(candidate_rank),
    "candidate_name": candidate_name,
    "status": status,
    "sequence_path": sequence_path,
    "structure_path": _to_project_rel(structure_path, project_root),
    "meta_path": _to_project_rel(meta_path, project_root),
    "duration_sec": duration_sec,
    "error_summary": error_summary if status == "failed" else "",
    "source_stage": source_stage,
    "source_name": source_name,
    "source_csv_path": _to_project_rel(source_csv_path, project_root) if source_csv_path else "",
    "source_row_idx": str(source_row_idx or ""),
    "source_sequence_field": source_sequence_field,
    "source_structure_field": source_structure_field,
    "native_score_name": native_score_name,
    "native_score_value": native_score_value,
    "notes": ";".join(note for note in notes if note),
  }


def _collect_mber_open(sample_id: str, native_root: Path, sample_output_dir: Path, project_root: Path, duration_sec: str) -> list[dict[str, str]]:
  accepted_csv = native_root / "accepted.csv"
  rows = _read_csv_rows(accepted_csv)
  if not rows:
    return []

  ordered_rows = _sort_rows(rows, rank_field_names=[], score_field_names=["i_ptm", "plddt", "ptm"])
  out: list[dict[str, str]] = []
  next_rank = 1
  for source_row_idx, row in ordered_rows:
    trajectory_name = (row.get("trajectory_name") or f"trajectory_{source_row_idx}").strip()
    binder_index = (row.get("binder_index") or str(source_row_idx)).strip()
    candidate_name = f"{trajectory_name}_binder-{binder_index}"
    sequence, sequence_field = _extract_sequence_from_row(row, preferred_fields=["binder_seq", "sequence"])

    structure_path = None
    structure_field = ""
    for field in ["relaxed_pdb_path", "complex_pdb_path"]:
      value = (row.get(field) or "").strip()
      if not value:
        continue
      structure_path = _resolve_candidate_path(value, accepted_csv.parent, [native_root, native_root / "Accepted"])
      if structure_path is not None:
        structure_field = field
        break
    if structure_path is None:
      structure_path = _glob_first(
        native_root,
        [
          f"Accepted/*{trajectory_name}_binder-{binder_index}_*relaxed.pdb",
          f"Accepted/*{trajectory_name}_binder-{binder_index}_*complex.pdb",
        ],
      )
      if structure_path is not None:
        structure_field = "glob_structure"

    score_name, score_value = _pick_native_score(row, ["i_ptm", "plddt", "ptm"])
    notes: list[str] = []
    if not sequence:
      notes.append("missing_sequence_in_csv")
    if structure_path is None:
      notes.append("missing_structure_path")

    out.append(
      _build_record(
        project_root=project_root,
        sample_output_dir=sample_output_dir,
        sample_id=sample_id,
        model="mber-open",
        candidate_rank=next_rank,
        candidate_name=candidate_name,
        sequence=sequence,
        structure_path=structure_path,
        duration_sec=duration_sec,
        source_stage="accepted",
        source_name=trajectory_name,
        source_csv_path=accepted_csv,
        source_row_idx=source_row_idx,
        source_sequence_field=sequence_field,
        source_structure_field=structure_field,
        native_score_name=score_name,
        native_score_value=score_value,
        error_summary="missing_candidate_artifacts" if not sequence and structure_path is None else "",
        notes=notes,
        raw_row=row,
      )
    )
    next_rank += 1
  return out


def _resolve_boltzgen_structure(native_root: Path, file_name: str) -> Path | None:
  file_name = (file_name or "").strip()
  if not file_name:
    return None
  patterns = [
    f"final_ranked_designs/**/before_refolding/*{file_name}",
    f"final_ranked_designs/**/*{file_name}",
    f"**/{file_name}",
  ]
  return _glob_first(native_root, patterns)


def _collect_boltzgen(sample_id: str, native_root: Path, sample_output_dir: Path, project_root: Path, duration_sec: str) -> list[dict[str, str]]:
  metrics_candidates = [native_root / "final_ranked_designs" / "all_designs_metrics.csv"]
  metrics_candidates.extend(sorted((native_root / "final_ranked_designs").glob("final_designs_metrics*.csv")))
  metrics_csv = next((path for path in metrics_candidates if path.exists()), None)
  if metrics_csv is None:
    return []
  rows = _read_csv_rows(metrics_csv)
  ordered_rows = _sort_rows(rows, rank_field_names=["final_rank", "Rank"], score_field_names=["quality_score", "design_to_target_iptm", "design_ptm"])
  out: list[dict[str, str]] = []
  for source_row_idx, row in ordered_rows:
    candidate_name = (row.get("id") or Path((row.get("file_name") or f"boltzgen_{source_row_idx}")).stem).strip()
    rank_text = (row.get("final_rank") or row.get("Rank") or "").strip()
    try:
      candidate_rank = int(float(rank_text)) if rank_text else len(out) + 1
    except ValueError:
      candidate_rank = len(out) + 1
    sequence, sequence_field = _extract_sequence_from_row(
      row,
      preferred_fields=["designed_chain_sequence", "designed_sequence", "Sequence", "sequence"],
    )
    structure_path = _resolve_boltzgen_structure(native_root, row.get("file_name", ""))
    score_name, score_value = _pick_native_score(row, ["quality_score", "design_to_target_iptm", "design_ptm", "pass_filters"])
    notes: list[str] = []
    if not sequence:
      notes.append("missing_sequence_in_metrics")
    if structure_path is None:
      notes.append("missing_structure_path")
    out.append(
      _build_record(
        project_root=project_root,
        sample_output_dir=sample_output_dir,
        sample_id=sample_id,
        model="boltzgen",
        candidate_rank=candidate_rank,
        candidate_name=candidate_name,
        sequence=sequence,
        structure_path=structure_path,
        duration_sec=duration_sec,
        source_stage="final_ranked_designs",
        source_name=candidate_name,
        source_csv_path=metrics_csv,
        source_row_idx=source_row_idx,
        source_sequence_field=sequence_field,
        source_structure_field="file_name",
        native_score_name=score_name,
        native_score_value=score_value,
        error_summary="missing_candidate_artifacts" if not sequence and structure_path is None else "",
        notes=notes,
        raw_row=row,
      )
    )
  return sorted(out, key=lambda row: (int(row["candidate_rank"]), row["candidate_name"]))


def _collect_germinal_rows_from_csv(
  *,
  rows: list[dict[str, str]],
  csv_path: Path,
  stage: str,
  native_root: Path,
  sample_id: str,
  sample_output_dir: Path,
  project_root: Path,
  duration_sec: str,
  seen_names: set[str],
  next_rank_start: int,
) -> tuple[list[dict[str, str]], int]:
  out: list[dict[str, str]] = []
  ordered_rows = _sort_rows(
    rows,
    rank_field_names=["rank"],
    score_field_names=["external_iptm", "i_ptm", "pdockq2", "plddt"],
  )
  next_rank = next_rank_start
  for source_row_idx, row in ordered_rows:
    candidate_name = (row.get("design_name") or row.get("Design") or f"germinal_{stage}_{source_row_idx}").strip()
    if candidate_name in seen_names:
      continue
    seen_names.add(candidate_name)
    sequence, sequence_field = _extract_sequence_from_row(row, preferred_fields=["sequence", "binder_sequence", "designed_sequence", "Sequence"])
    structure_path = None
    structure_field = ""
    for field in ["final_structure_path", "structure_path", "pdb_path", "structure", "pdb"]:
      value = (row.get(field) or "").strip()
      if not value:
        continue
      structure_path = _resolve_candidate_path(value, csv_path.parent, [native_root, csv_path.parent.parent])
      if structure_path is not None:
        structure_field = field
        break
    if structure_path is None:
      structure_path = _glob_first(native_root, [f"**/structures/*{candidate_name}*.pdb", f"**/structures/*{candidate_name}*.cif"])
      if structure_path is not None:
        structure_field = "glob_structure"
    score_name, score_value = _pick_native_score(row, ["external_iptm", "i_ptm", "pdockq2", "plddt"])
    notes: list[str] = []
    if not sequence:
      notes.append("missing_sequence_in_csv")
    if structure_path is None:
      notes.append("missing_structure_path")
    out.append(
      _build_record(
        project_root=project_root,
        sample_output_dir=sample_output_dir,
        sample_id=sample_id,
        model="germinal",
        candidate_rank=next_rank,
        candidate_name=candidate_name,
        sequence=sequence,
        structure_path=structure_path,
        duration_sec=duration_sec,
        source_stage=stage,
        source_name=candidate_name,
        source_csv_path=csv_path,
        source_row_idx=source_row_idx,
        source_sequence_field=sequence_field,
        source_structure_field=structure_field,
        native_score_name=score_name,
        native_score_value=score_value,
        error_summary="missing_candidate_artifacts" if not sequence and structure_path is None else "",
        notes=notes,
        raw_row=row,
      )
    )
    next_rank += 1
  return out, next_rank


def _collect_germinal_structure_only(
  *,
  native_root: Path,
  sample_id: str,
  sample_output_dir: Path,
  project_root: Path,
  duration_sec: str,
  seen_names: set[str],
  next_rank_start: int,
) -> tuple[list[dict[str, str]], int]:
  out: list[dict[str, str]] = []
  next_rank = next_rank_start
  for stage in ["accepted", "redesign_candidates", "trajectories"]:
    structure_paths = sorted(native_root.glob(f"**/{stage}/structures/*.pdb"))
    structure_paths += sorted(native_root.glob(f"**/{stage}/structures/*.cif"))
    for structure_path in structure_paths:
      candidate_name = structure_path.stem.strip() or f"germinal_{stage}_{next_rank}"
      if candidate_name in seen_names:
        continue
      seen_names.add(candidate_name)
      out.append(
        _build_record(
          project_root=project_root,
          sample_output_dir=sample_output_dir,
          sample_id=sample_id,
          model="germinal",
          candidate_rank=next_rank,
          candidate_name=candidate_name,
          sequence="",
          structure_path=structure_path,
          duration_sec=duration_sec,
          source_stage=stage,
          source_name=candidate_name,
          source_csv_path=None,
          source_row_idx=None,
          source_sequence_field="",
          source_structure_field="glob_structure_only",
          native_score_name="",
          native_score_value="",
          error_summary="",
          notes=["structure_only_fallback", "missing_sequence_in_csv"],
          raw_row=None,
        )
      )
      next_rank += 1
  return out, next_rank


def _collect_germinal(sample_id: str, native_root: Path, sample_output_dir: Path, project_root: Path, duration_sec: str) -> list[dict[str, str]]:
  all_csvs = sorted(native_root.glob("**/all_trajectories.csv"))
  seen_names: set[str] = set()
  next_rank = 1
  out: list[dict[str, str]] = []
  if all_csvs:
    for csv_path in all_csvs:
      rows, next_rank = _collect_germinal_rows_from_csv(
        rows=_read_csv_rows(csv_path),
        csv_path=csv_path,
        stage="all_trajectories",
        native_root=native_root,
        sample_id=sample_id,
        sample_output_dir=sample_output_dir,
        project_root=project_root,
        duration_sec=duration_sec,
        seen_names=seen_names,
        next_rank_start=next_rank,
      )
      out.extend(rows)

  stage_specs = [
    ("accepted", sorted(native_root.glob("**/accepted/designs.csv"))),
    ("redesign_candidates", sorted(native_root.glob("**/redesign_candidates/designs.csv"))),
    ("trajectories", sorted(native_root.glob("**/trajectories/designs.csv"))),
  ]
  for stage, csv_paths in stage_specs:
    for csv_path in csv_paths:
      rows, next_rank = _collect_germinal_rows_from_csv(
        rows=_read_csv_rows(csv_path),
        csv_path=csv_path,
        stage=stage,
        native_root=native_root,
        sample_id=sample_id,
        sample_output_dir=sample_output_dir,
        project_root=project_root,
        duration_sec=duration_sec,
        seen_names=seen_names,
        next_rank_start=next_rank,
      )
      out.extend(rows)

  fallback_rows, next_rank = _collect_germinal_structure_only(
    native_root=native_root,
    sample_id=sample_id,
    sample_output_dir=sample_output_dir,
    project_root=project_root,
    duration_sec=duration_sec,
    seen_names=seen_names,
    next_rank_start=next_rank,
  )
  out.extend(fallback_rows)
  return out


def _split_globs(raw: str, default_patterns: list[str]) -> list[str]:
  text = (raw or "").strip()
  if not text:
    return default_patterns
  return [part.strip() for part in text.split(":") if part.strip()]


def _match_sequence_file(structure_path: Path, sequence_files: list[Path]) -> Path | None:
  stem = structure_path.stem
  exact = [path for path in sequence_files if path.stem == stem]
  if exact:
    return exact[0]
  contains = [path for path in sequence_files if stem in path.stem or path.stem in stem]
  if contains:
    return contains[0]
  return None


def _rank_from_rfantibody_name(name: str, fallback: int) -> int:
  match = re.search(r"samples_design_(\d+)_dldesign_(\d+)_best", name)
  if match:
    try:
      primary = int(match.group(1))
      secondary = int(match.group(2))
      return primary * 100 + secondary + 1
    except ValueError:
      return fallback
  return fallback


def _collect_rfantibody(sample_id: str, native_root: Path, sample_output_dir: Path, project_root: Path, duration_sec: str) -> list[dict[str, str]]:
  structure_patterns = _split_globs(
    os.environ.get("RFANTIBODY_CANDIDATE_GLOBS", ""),
    [
      "final_designs/*.pdb",
      "selected_designs/*.pdb",
      "predictions/*.pdb",
      "native_run*.pdb",
      "native_run*.cif",
      "**/*best.pdb",
      "**/*.pdb",
      "**/*.cif",
    ],
  )
  sequence_patterns = _split_globs(
    os.environ.get("RFANTIBODY_SEQUENCE_GLOBS", ""),
    ["native_run*.fa", "native_run*.faa", "native_run*.fasta", "native_run*.fas", "native_run*.seq", "**/*.fa", "**/*.faa", "**/*.fasta", "**/*.fas", "**/*.seq"],
  )

  search_roots = [native_root]
  if sample_output_dir != native_root:
    search_roots.append(sample_output_dir)

  structure_paths: list[Path] = []
  sequence_files: list[Path] = []
  seen_structures: set[str] = set()
  seen_sequences: set[str] = set()
  for root in search_roots:
    for path in _glob_all(root, structure_patterns):
      if path.suffix.lower() not in STRUCTURE_EXTENSIONS:
        continue
      key = str(path.resolve()) if path.exists() else str(path)
      if key in seen_structures:
        continue
      seen_structures.add(key)
      structure_paths.append(path)
    for path in _glob_all(root, sequence_patterns):
      if path.suffix.lower() not in FASTA_EXTENSIONS:
        continue
      key = str(path.resolve()) if path.exists() else str(path)
      if key in seen_sequences:
        continue
      seen_sequences.add(key)
      sequence_files.append(path)

  out: list[dict[str, str]] = []
  seen_names: set[str] = set()
  for idx, structure_path in enumerate(structure_paths, start=1):
    candidate_name = structure_path.stem
    if candidate_name in seen_names:
      continue
    seen_names.add(candidate_name)
    sequence_file = _match_sequence_file(structure_path, sequence_files)
    sequence = ""
    sequence_field = ""
    if sequence_file is not None:
      lines = sequence_file.read_text(encoding="utf-8").splitlines()
      sequence = _normalize_sequence("".join(line.strip() for line in lines if not line.startswith(">")))
      sequence_field = "existing_sequence_file"
    notes: list[str] = []
    if not sequence:
      notes.append("missing_explicit_sequence_file")
    candidate_rank = _rank_from_rfantibody_name(candidate_name, idx)
    out.append(
      _build_record(
        project_root=project_root,
        sample_output_dir=sample_output_dir,
        sample_id=sample_id,
        model="RFantibody",
        candidate_rank=candidate_rank,
        candidate_name=candidate_name,
        sequence=sequence,
        structure_path=structure_path,
        duration_sec=duration_sec,
        source_stage="native_outputs",
        source_name=candidate_name,
        source_csv_path=None,
        source_row_idx=None,
        source_sequence_field=sequence_field,
        source_structure_field="glob_structure",
        native_score_name="",
        native_score_value="",
        error_summary="",
        notes=notes,
        raw_row={
          "structure_path": str(structure_path),
          "sequence_file": str(sequence_file) if sequence_file else "",
        },
      )
    )
  return sorted(out, key=lambda row: (int(row["candidate_rank"]), row["candidate_name"]))


def collect_candidates(
  *,
  model: str,
  sample_id: str,
  sample_input_dir: Path | None,
  sample_output_dir: Path,
  project_root: Path,
  runner_exit_code: int,
  duration_sec: str,
) -> tuple[list[dict[str, str]], dict]:
  native_root = sample_output_dir / "native_run"
  records: list[dict[str, str]]
  model_key = model.strip()
  if model_key == "boltzgen":
    records = _collect_boltzgen(sample_id, native_root, sample_output_dir, project_root, duration_sec)
  elif model_key == "germinal":
    records = _collect_germinal(sample_id, native_root, sample_output_dir, project_root, duration_sec)
  elif model_key == "mber-open":
    records = _collect_mber_open(sample_id, native_root, sample_output_dir, project_root, duration_sec)
  elif model_key == "RFantibody":
    records = _collect_rfantibody(sample_id, native_root, sample_output_dir, project_root, duration_sec)
  else:
    records = []

  notes: list[str] = []
  sample_status = "ok"
  error_summary = ""

  if not records:
    sample_status = "failed"
    if runner_exit_code != 0:
      notes.append(f"runner_exit_nonzero:{runner_exit_code}")
      error_summary = f"runner_exit_nonzero:{runner_exit_code}"
    else:
      notes.append("no_candidates_found")
      error_summary = "no_candidates_found"
    failure_record = _build_record(
      project_root=project_root,
      sample_output_dir=sample_output_dir,
      sample_id=sample_id,
      model=model_key,
      candidate_rank=1,
      candidate_name="run_failure",
      sequence="",
      structure_path=None,
      duration_sec=duration_sec,
      source_stage="run_failure",
      source_name="run_failure",
      source_csv_path=None,
      source_row_idx=None,
      source_sequence_field="",
      source_structure_field="",
      native_score_name="",
      native_score_value="",
      error_summary=f"runner_exit_nonzero:{runner_exit_code}" if runner_exit_code != 0 else "no_candidates_found",
      notes=[f"runner_exit_nonzero:{runner_exit_code}" if runner_exit_code != 0 else "no_candidates_found"],
      raw_row={},
    )
    records = [failure_record]
  elif runner_exit_code != 0:
    notes.append(f"runner_exit_nonzero:{runner_exit_code}")
    sample_status = "partial"
    error_summary = f"runner_exit_nonzero:{runner_exit_code}"

  candidate_manifest = sample_output_dir / "candidate_manifest.csv"
  candidate_manifest.parent.mkdir(parents=True, exist_ok=True)
  with candidate_manifest.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
    writer.writeheader()
    writer.writerows(records)

  run_meta = {
    "status": sample_status,
    "model": model_key,
    "sample_id": sample_id,
    "runner_exit_code": runner_exit_code,
    "duration_sec": duration_sec,
    "n_candidates": len([row for row in records if row["status"] == "ok"]),
    "candidate_manifest_path": _to_project_rel(candidate_manifest, project_root),
    "error_summary": error_summary,
    "notes": notes,
    "sample_input_dir": _to_project_rel(sample_input_dir, project_root) if sample_input_dir else "",
    "sample_output_dir": _to_project_rel(sample_output_dir, project_root),
  }
  _write_json(sample_output_dir / "run_meta.json", run_meta)
  return records, run_meta


def main() -> int:
  parser = argparse.ArgumentParser(description="Collect all generated candidates into a standardized manifest.")
  parser.add_argument("--print-fieldnames", action="store_true", help="Print manifest CSV header and exit.")
  parser.add_argument("--model", type=str, default="")
  parser.add_argument("--sample-id", type=str, default="")
  parser.add_argument("--sample-input-dir", type=Path, default=None)
  parser.add_argument("--sample-output-dir", type=Path, default=None)
  parser.add_argument("--project-root", type=Path, default=Path.cwd())
  parser.add_argument("--runner-exit-code", type=int, default=0)
  parser.add_argument("--duration-sec", type=str, default="")
  args = parser.parse_args()

  if args.print_fieldnames:
    return _print_fieldnames()

  if not args.model or not args.sample_id or args.sample_output_dir is None:
    parser.error("--model、--sample-id、--sample-output-dir 为必填参数")

  collect_candidates(
    model=args.model,
    sample_id=args.sample_id,
    sample_input_dir=args.sample_input_dir,
    sample_output_dir=args.sample_output_dir,
    project_root=args.project_root,
    runner_exit_code=args.runner_exit_code,
    duration_sec=args.duration_sec,
  )
  return 0


if __name__ == "__main__":
  raise SystemExit(main())