#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
from pathlib import Path

from evaluation.schema import EXTERNAL_CANDIDATE_FIELDS, EXTERNAL_SEED_FIELDS


def _load_csv(path: Path) -> list[dict[str, str]]:
  if not path.exists():
    return []
  with path.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)


def _safe_read_json(path: Path) -> dict:
  if not path.exists():
    return {}
  try:
    return json.loads(path.read_text(encoding="utf-8"))
  except Exception:
    return {}


def _resolve_path(repo_root: Path, raw_path: str, *, base_dir: Path | None = None) -> str:
  text = (raw_path or "").strip()
  if not text:
    return ""
  path = Path(text)
  if not path.is_absolute():
    path = (base_dir / path) if base_dir is not None else (repo_root / path)
  return str(path.resolve())


def _to_int(value: str) -> int | None:
  text = (value or "").strip()
  if not text:
    return None
  try:
    return int(float(text))
  except ValueError:
    return None


def _result_value(payload: dict, *keys: str) -> str:
  for key in keys:
    value = payload.get(key)
    if value is None:
      continue
    text = str(value).strip()
    if text:
      return text
  return ""


def _manifest_by_candidate(judge_manifest: Path | None) -> dict[str, list[dict[str, str]]]:
  out: dict[str, list[dict[str, str]]] = {}
  if judge_manifest is None or not judge_manifest.exists():
    return out
  for row in _load_csv(judge_manifest):
    candidate_id = (row.get("candidate_id") or "").strip()
    if candidate_id:
      out.setdefault(candidate_id, []).append(row)
  return out


def _base_seed_row(candidate: dict[str, str]) -> dict[str, str]:
  row = {k: "" for k in EXTERNAL_SEED_FIELDS}
  for key in EXTERNAL_CANDIDATE_FIELDS:
    if key in row:
      row[key] = candidate.get(key, "") or ""
  for key in [
    "candidate_id",
    "run_id",
    "run_tag",
    "model",
    "manifest_row_idx",
    "source_candidate_rank_in_model",
    "candidate_rank_in_model",
    "sample_id",
    "pred_sequence",
    "pred_seq_len",
    "aa_valid_ratio",
    "cys_count",
    "reference_complex_path",
    "target_structure_path",
    "target_chain_id",
    "designed_antibody_chain_id",
    "epitope_path",
    "hotspot_count",
    "hotspot_string",
    "source_duration_sec",
    "num_missing_residues",
  ]:
    row[key] = candidate.get(key, "") or ""
  return row


def _infer_structure_path(output_dir: Path) -> Path | None:
  for pattern in ("*.cif", "*.mmcif", "*.pdb"):
    matches = sorted(path for path in output_dir.rglob(pattern) if path.is_file())
    if matches:
      return matches[0]
  return None


def _infer_confidence_path(output_dir: Path, result_json: Path) -> Path | None:
  excluded = {result_json.resolve()}
  for pattern in ("*.json", "*.npz", "*.pkl", "*.pt", "*.csv"):
    matches = sorted(
      path
      for path in output_dir.rglob(pattern)
      if path.is_file() and path.resolve() not in excluded
    )
    if matches:
      return matches[0]
  return None


def _render_command(template: str, replacements: dict[str, str]) -> str:
  rendered = template
  for placeholder, value in replacements.items():
    rendered = rendered.replace(placeholder, shlex.quote(value))
  return rendered


def _build_imported_row(
  *,
  candidate: dict[str, str],
  imported: dict[str, str],
  repo_root: Path,
  base_dir: Path | None,
) -> dict[str, str]:
  row = _base_seed_row(candidate)
  row["judge_name"] = (imported.get("judge_name") or candidate.get("judge_name") or "chai1_primary").strip()
  row["judge_seed"] = (imported.get("judge_seed") or "0").strip()
  row["judge_mode"] = "imported_manifest"
  row["external_status"] = (imported.get("external_status") or imported.get("status") or "").strip()
  row["external_error_summary"] = (imported.get("external_error_summary") or imported.get("error_summary") or "").strip()
  row["external_structure_path"] = _resolve_path(
    repo_root,
    imported.get("external_structure_path") or imported.get("structure_path") or "",
    base_dir=base_dir,
  )
  row["external_confidence_path"] = _resolve_path(
    repo_root,
    imported.get("external_confidence_path") or imported.get("confidence_path") or "",
    base_dir=base_dir,
  )
  row["external_iptm"] = _result_value(imported, "external_iptm", "iptm")
  row["external_ipae"] = _result_value(imported, "external_ipae", "ipae")
  row["external_plddt_binder"] = _result_value(imported, "external_plddt_binder", "plddt_binder")
  if (imported.get("target_chain_id") or "").strip():
    row["target_chain_id"] = (imported.get("target_chain_id") or "").strip()
  if (imported.get("designed_antibody_chain_id") or "").strip():
    row["designed_antibody_chain_id"] = (imported.get("designed_antibody_chain_id") or "").strip()
  elif (imported.get("antibody_chain_id") or "").strip():
    row["designed_antibody_chain_id"] = (imported.get("antibody_chain_id") or "").strip()
  return row


def _build_execution_row(
  *,
  candidate: dict[str, str],
  judge_seed: int,
  judge_name: str,
  repo_root: Path,
  output_dir: Path,
  result_json: Path,
  returncode: int,
  forced_error: str = "",
) -> dict[str, str]:
  payload = _safe_read_json(result_json)
  row = _base_seed_row(candidate)
  row["judge_name"] = judge_name
  row["judge_seed"] = str(judge_seed)
  row["judge_mode"] = "executed_command"

  structure_path = _resolve_path(
    repo_root,
    _result_value(payload, "external_structure_path", "structure_path"),
    base_dir=output_dir,
  )
  if not structure_path:
    inferred_structure = _infer_structure_path(output_dir)
    structure_path = str(inferred_structure.resolve()) if inferred_structure else ""

  confidence_path = _resolve_path(
    repo_root,
    _result_value(payload, "external_confidence_path", "confidence_path"),
    base_dir=output_dir,
  )
  if not confidence_path:
    inferred_confidence = _infer_confidence_path(output_dir, result_json)
    confidence_path = str(inferred_confidence.resolve()) if inferred_confidence else ""

  status = _result_value(payload, "external_status", "status")
  if status == "ok" and not structure_path:
    status = "missing_external_structure"
  if not status:
    if structure_path:
      status = "ok"
    elif returncode == 0:
      status = "missing_external_structure"
    else:
      status = "failed"

  error_summary = forced_error or _result_value(payload, "external_error_summary", "error_summary")
  if not error_summary and status != "ok":
    error_summary = f"judge_cmd_failed:{returncode}" if returncode != 0 else "missing_external_structure"

  row["external_status"] = status
  row["external_error_summary"] = error_summary
  row["external_structure_path"] = structure_path
  row["external_confidence_path"] = confidence_path
  row["external_iptm"] = _result_value(payload, "external_iptm", "iptm")
  row["external_ipae"] = _result_value(payload, "external_ipae", "ipae")
  row["external_plddt_binder"] = _result_value(payload, "external_plddt_binder", "plddt_binder")
  if _result_value(payload, "target_chain_id"):
    row["target_chain_id"] = _result_value(payload, "target_chain_id")
  if _result_value(payload, "designed_antibody_chain_id"):
    row["designed_antibody_chain_id"] = _result_value(payload, "designed_antibody_chain_id")
  elif _result_value(payload, "antibody_chain_id"):
    row["designed_antibody_chain_id"] = _result_value(payload, "antibody_chain_id")
  return row


def _run_judge_command(
  *,
  candidate: dict[str, str],
  judge_seed: int,
  judge_name: str,
  judge_cmd: str,
  repo_root: Path,
  output_dir: Path,
) -> dict[str, str]:
  output_dir.mkdir(parents=True, exist_ok=True)
  stdout_log = output_dir / "stdout.log"
  stderr_log = output_dir / "stderr.log"
  result_json = output_dir / "result.json"

  replacements = {
    "__CANDIDATE_ID__": candidate.get("candidate_id", "") or "",
    "__SAMPLE_ID__": candidate.get("sample_id", "") or "",
    "__MODEL__": candidate.get("model", "") or "",
    "__PRED_SEQUENCE__": candidate.get("pred_sequence", "") or "",
    "__SOURCE_SEQUENCE_PATH__": candidate.get("source_sequence_path", "") or "",
    "__SOURCE_STRUCTURE_PATH__": candidate.get("source_structure_path", "") or "",
    "__TARGET_STRUCTURE_PATH__": candidate.get("target_structure_path", "") or "",
    "__TARGET_CHAIN_ID__": candidate.get("target_chain_id", "") or "",
    "__DESIGNED_ANTIBODY_CHAIN_ID__": candidate.get("designed_antibody_chain_id", "") or "",
    "__HOTSPOT_STRING__": candidate.get("hotspot_string", "") or "",
    "__JUDGE_SEED__": str(judge_seed),
    "__JUDGE_NAME__": judge_name,
    "__JOB_JSON__": candidate.get("judge_input_path", "") or "",
    "__OUTPUT_DIR__": str(output_dir.resolve()),
    "__RESULT_JSON__": str(result_json.resolve()),
  }
  command = _render_command(judge_cmd, replacements)
  env = os.environ.copy()
  env.update(
    {
      "PRIMARY_JUDGE_CANDIDATE_ID": replacements["__CANDIDATE_ID__"],
      "PRIMARY_JUDGE_SAMPLE_ID": replacements["__SAMPLE_ID__"],
      "PRIMARY_JUDGE_MODEL": replacements["__MODEL__"],
      "PRIMARY_JUDGE_PRED_SEQUENCE": replacements["__PRED_SEQUENCE__"],
      "PRIMARY_JUDGE_SOURCE_SEQUENCE_PATH": replacements["__SOURCE_SEQUENCE_PATH__"],
      "PRIMARY_JUDGE_SOURCE_STRUCTURE_PATH": replacements["__SOURCE_STRUCTURE_PATH__"],
      "PRIMARY_JUDGE_TARGET_STRUCTURE_PATH": replacements["__TARGET_STRUCTURE_PATH__"],
      "PRIMARY_JUDGE_TARGET_CHAIN_ID": replacements["__TARGET_CHAIN_ID__"],
      "PRIMARY_JUDGE_DESIGNED_ANTIBODY_CHAIN_ID": replacements["__DESIGNED_ANTIBODY_CHAIN_ID__"],
      "PRIMARY_JUDGE_HOTSPOT_STRING": replacements["__HOTSPOT_STRING__"],
      "PRIMARY_JUDGE_SEED": replacements["__JUDGE_SEED__"],
      "PRIMARY_JUDGE_NAME": replacements["__JUDGE_NAME__"],
      "PRIMARY_JUDGE_JOB_JSON": replacements["__JOB_JSON__"],
      "PRIMARY_JUDGE_OUTPUT_DIR": replacements["__OUTPUT_DIR__"],
      "PRIMARY_JUDGE_RESULT_JSON": replacements["__RESULT_JSON__"],
    }
  )

  try:
    completed = subprocess.run(
      command,
      shell=True,
      executable="/bin/bash",
      cwd=str(repo_root),
      env=env,
      capture_output=True,
      text=True,
      check=False,
    )
    stdout_log.write_text(completed.stdout or "", encoding="utf-8")
    stderr_log.write_text(completed.stderr or "", encoding="utf-8")
    return _build_execution_row(
      candidate=candidate,
      judge_seed=judge_seed,
      judge_name=judge_name,
      repo_root=repo_root,
      output_dir=output_dir,
      result_json=result_json,
      returncode=completed.returncode,
    )
  except Exception as exc:
    stdout_log.write_text("", encoding="utf-8")
    stderr_log.write_text(f"{type(exc).__name__}: {exc}\n", encoding="utf-8")
    return _build_execution_row(
      candidate=candidate,
      judge_seed=judge_seed,
      judge_name=judge_name,
      repo_root=repo_root,
      output_dir=output_dir,
      result_json=result_json,
      returncode=1,
      forced_error=f"judge_cmd_exception:{type(exc).__name__}:{exc}",
    )


def run_primary_judge(
  *,
  candidates_csv: Path,
  out_csv: Path,
  judge_manifest: Path | None,
  fallback_to_source_structure: bool,
  judge_cmd: str | None = None,
  run_dir: Path | None = None,
) -> Path:
  repo_root = Path(__file__).resolve().parents[2]
  judge_cmd = (judge_cmd or os.environ.get("PRIMARY_JUDGE_CMD") or "").strip() or None
  run_dir = run_dir or (out_csv.parent / "primary_judge_runs")
  candidates = _load_csv(candidates_csv)
  judge_rows_by_candidate = _manifest_by_candidate(judge_manifest)
  imported_base_dir = judge_manifest.parent if judge_manifest else None
  out_rows: list[dict[str, str]] = []

  for candidate in candidates:
    candidate_id = (candidate.get("candidate_id") or "").strip()
    imported_rows = judge_rows_by_candidate.get(candidate_id, [])

    if imported_rows:
      for imported in imported_rows:
        out_rows.append(
          _build_imported_row(
            candidate=candidate,
            imported=imported,
            repo_root=repo_root,
            base_dir=imported_base_dir,
          )
        )
      continue

    if judge_cmd:
      n_requested_seeds = max(1, _to_int(candidate.get("n_requested_seeds", "")) or 1)
      candidate_slug = candidate_id.replace(":", "__").replace("/", "_")
      for judge_seed in range(n_requested_seeds):
        out_rows.append(
          _run_judge_command(
            candidate=candidate,
            judge_seed=judge_seed,
            judge_name=(candidate.get("judge_name") or "chai1_primary").strip(),
            judge_cmd=judge_cmd,
            repo_root=repo_root,
            output_dir=run_dir / candidate_slug / f"seed_{judge_seed:03d}",
          )
        )
      continue

    row = _base_seed_row(candidate)
    if fallback_to_source_structure:
      row["judge_name"] = (candidate.get("judge_name") or "chai1_primary").strip()
      row["judge_seed"] = "0"
      row["judge_mode"] = "source_structure_fallback"
      row["external_structure_path"] = _resolve_path(repo_root, candidate.get("source_structure_path") or "")
      row["external_status"] = "ok" if row["external_structure_path"] else "failed"
      row["external_error_summary"] = "" if row["external_status"] == "ok" else "missing_source_structure"
    else:
      row["judge_name"] = (candidate.get("judge_name") or "chai1_primary").strip()
      row["judge_seed"] = "0"
      row["judge_mode"] = "pending_import"
      row["external_structure_path"] = ""
      row["external_status"] = "missing_judge_output"
      row["external_error_summary"] = "missing_judge_output"
    out_rows.append(row)

  _write_csv(out_csv, EXTERNAL_SEED_FIELDS, out_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="导入主 judge 结果，或用已有结构作为 Phase1 占位输入")
  parser.add_argument(
    "--candidates-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_candidates.csv"),
    help="候选表路径",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="seed 级主表（含 judge 输出）",
  )
  parser.add_argument(
    "--judge-manifest",
    type=Path,
    default=None,
    help="外部 judge 结果 manifest，可选字段 candidate_id/judge_seed/external_structure_path/status",
  )
  parser.add_argument(
    "--fallback-to-source-structure",
    action="store_true",
    help="若没有 judge manifest，则回退使用模型原始结构作为 Phase1 占位输入",
  )
  parser.add_argument(
    "--judge-cmd",
    type=str,
    default=None,
    help="直接执行主 judge 的命令模板，可用占位符如 __JOB_JSON__/__JUDGE_SEED__/__OUTPUT_DIR__/__RESULT_JSON__",
  )
  parser.add_argument(
    "--run-dir",
    type=Path,
    default=None,
    help="主 judge 执行目录，默认 out-csv 同级 primary_judge_runs",
  )
  args = parser.parse_args()

  out_csv = run_primary_judge(
    candidates_csv=args.candidates_csv,
    out_csv=args.out_csv,
    judge_manifest=args.judge_manifest,
    fallback_to_source_structure=args.fallback_to_source_structure,
    judge_cmd=args.judge_cmd,
    run_dir=args.run_dir,
  )
  print(f"[OK] 主 judge seed 表: {out_csv}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())