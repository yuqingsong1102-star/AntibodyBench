#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path


READY_FIELDS = [
  "sample_id",
  "pdb_id",
  "reference_structure_path",
  "design_structure_path",
  "design_structure_source",
  "antigen_chain",
  "antibody_chain",
  "reference_complex_path",
  "reference_complex_status",
  "reference_complex_note",
]


def _safe(text: str) -> str:
  return (text or "").strip()


def _base_pdb_id(raw: str) -> str:
  return _safe(raw).split("-")[0].lower()


def _is_nonempty_file(path: Path) -> bool:
  return path.exists() and path.is_file() and path.stat().st_size > 20


def _download(url: str, out_path: Path, timeout_sec: int) -> bool:
  try:
    with urllib.request.urlopen(url, timeout=timeout_sec) as resp:
      data = resp.read()
  except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
    return False
  if not data:
    return False
  out_path.write_bytes(data)
  return True


def _run_python_step(script_path: Path, args: list[str]) -> None:
  cmd = [sys.executable, str(script_path), *args]
  print(f"[INFO] step={script_path.name}")
  try:
    subprocess.run(cmd, check=True)
  except subprocess.CalledProcessError as exc:
    raise SystemExit(f"[ERROR] step failed: {script_path.name} (exit={exc.returncode})") from exc


def _load_rows(path: Path) -> list[dict[str, str]]:
  with path.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _build_base_rows(raw_index: Path, reference_pdb_dir: Path) -> list[dict[str, str]]:
  base_rows: list[dict[str, str]] = []
  for raw_row in _load_rows(raw_index):
    raw_pdb_id = _safe(raw_row.get("pdb_id", "") or raw_row.get("PDB ID", ""))
    antibody_chain = _safe(raw_row.get("antibody_chain", ""))
    antigen_chain = _safe(raw_row.get("antigen_chain", ""))
    if not raw_pdb_id or not antibody_chain or not antigen_chain:
      continue

    pdb_id = _base_pdb_id(raw_pdb_id)
    sample_id = f"{pdb_id}_{antibody_chain}_{antigen_chain}"
    base_rows.append(
      {
        "sample_id": sample_id,
        "pdb_id": pdb_id,
        "reference_structure_path": str((reference_pdb_dir / f"{pdb_id}.pdb").resolve()),
        "design_structure_path": str((reference_pdb_dir / f"{pdb_id}.pdb").resolve()),
        "design_structure_source": "reference_structure_path",
        "antigen_chain": antigen_chain,
        "antibody_chain": antibody_chain,
        "reference_complex_path": "",
        "reference_complex_status": "",
        "reference_complex_note": "",
      }
    )
  return base_rows


def _fill_reference_complex_paths(rows: list[dict[str, str]], complex_dir: Path, timeout_sec: int) -> tuple[int, int]:
  complex_dir.mkdir(parents=True, exist_ok=True)
  cache: dict[str, tuple[str, str, str]] = {}
  ok = 0
  fail = 0

  for row in rows:
    pdb_id = _safe(row.get("pdb_id", ""))
    if not pdb_id:
      row["reference_complex_status"] = "fail"
      row["reference_complex_note"] = "missing_pdb_id"
      fail += 1
      continue

    if pdb_id in cache:
      path_value, status_value, note_value = cache[pdb_id]
      row["reference_complex_path"] = path_value
      row["reference_complex_status"] = status_value
      row["reference_complex_note"] = note_value
      if status_value == "ok":
        ok += 1
      else:
        fail += 1
      continue

    out_pdb = complex_dir / f"{pdb_id}.pdb"
    out_cif = complex_dir / f"{pdb_id}.cif"
    if _is_nonempty_file(out_pdb):
      cache[pdb_id] = (str(out_pdb.resolve()), "ok", "")
    elif _is_nonempty_file(out_cif):
      cache[pdb_id] = (str(out_cif.resolve()), "ok", "downloaded_cif")
    else:
      pdb_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
      cif_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
      if _download(pdb_url, out_pdb, timeout_sec) and _is_nonempty_file(out_pdb):
        cache[pdb_id] = (str(out_pdb.resolve()), "ok", "")
      elif _download(cif_url, out_cif, timeout_sec) and _is_nonempty_file(out_cif):
        cache[pdb_id] = (str(out_cif.resolve()), "ok", "downloaded_cif")
      else:
        cache[pdb_id] = ("", "fail", "download_failed")

    path_value, status_value, note_value = cache[pdb_id]
    row["reference_complex_path"] = path_value
    row["reference_complex_status"] = status_value
    row["reference_complex_note"] = note_value
    if status_value == "ok":
      ok += 1
    else:
      fail += 1

  return ok, fail


def _write_ready_rows(path: Path, rows: list[dict[str, str]]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=READY_FIELDS)
    writer.writeheader()
    for row in rows:
      writer.writerow({key: _safe(str(row.get(key, ""))) for key in READY_FIELDS})


def _select_ready_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
  ready_rows: list[dict[str, str]] = []
  seen: set[str] = set()
  for row in rows:
    sample_id = _safe(row.get("sample_id", ""))
    if not sample_id or sample_id in seen:
      continue
    if _safe(row.get("reference_complex_status", "")) != "ok":
      continue
    seen.add(sample_id)
    ready_rows.append(row)
  return ready_rows


def _keep_rows_with_target_paths(rows: list[dict[str, str]]) -> list[dict[str, str]]:
  out: list[dict[str, str]] = []
  for row in rows:
    target_path = Path(_safe(row.get("design_structure_path", "")) or _safe(row.get("reference_structure_path", "")))
    if _is_nonempty_file(target_path):
      out.append(row)
  return out


def main() -> int:
  parser = argparse.ArgumentParser(description="Prepare the four model input trees in one command.")
  parser.add_argument(
    "--raw-index",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "dataset_index.csv",
    help="Source dataset index CSV.",
  )
  parser.add_argument(
    "--reference-pdb-dir",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "pdbs",
    help="Directory containing target/reference PDB files (auto-downloaded if missing).",
  )
  parser.add_argument(
    "--complex-dir",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "complexes",
    help="Directory containing downloaded full complexes.",
  )
  parser.add_argument(
    "--ready-csv",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "dataset_index_ready.csv",
    help="Final ready index consumed by downstream build/run/eval scripts.",
  )
  parser.add_argument(
    "--epitope-dir",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "epitopes",
    help="Output directory for per-sample epitope JSON files.",
  )
  parser.add_argument(
    "--native-root",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "native_inputs",
    help="Output directory for native model inputs.",
  )
  parser.add_argument(
    "--models-root",
    type=Path,
    default=Path(__file__).resolve().parents[3] / "models",
    help="Root directory containing the four upstream model repos.",
  )
  parser.add_argument(
    "--timeout-sec",
    type=int,
    default=8,
    help="HTTP timeout for downloading missing full complexes.",
  )
  parser.add_argument(
    "--cutoff",
    type=float,
    default=5.0,
    help="Heavy-atom distance cutoff for epitope extraction.",
  )
  parser.add_argument(
    "--no-repair-target-pdbs",
    action="store_true",
    help="Disable target-chain repair from full complexes during epitope extraction.",
  )
  parser.add_argument(
    "--no-crop-targets",
    action="store_true",
    help="Disable automatic cropped design-target generation for large antigens.",
  )
  parser.add_argument(
    "--cropped-target-dir",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "cropped_targets",
    help="Output directory for automatically cropped design targets.",
  )
  parser.add_argument(
    "--crop-radius",
    type=float,
    default=18.0,
    help="3D radius in Angstrom around hotspot residues to retain when cropping large targets.",
  )
  parser.add_argument(
    "--crop-padding-residues",
    type=int,
    default=8,
    help="Sequence padding to retain on each side of selected hotspot-neighbor residues.",
  )
  parser.add_argument(
    "--crop-min-residues",
    type=int,
    default=350,
    help="Only auto-crop targets whose antigen residue count exceeds this threshold.",
  )
  args = parser.parse_args()

  repo_root = Path(__file__).resolve().parents[2]
  data_prep_dir = Path(__file__).resolve().parent
  if not args.raw_index.exists():
    raise SystemExit(f"[ERROR] raw index not found: {args.raw_index}")

  args.reference_pdb_dir.mkdir(parents=True, exist_ok=True)
  rows = _build_base_rows(args.raw_index, args.reference_pdb_dir)
  ok_complexes, fail_complexes = _fill_reference_complex_paths(rows, args.complex_dir, args.timeout_sec)

  ready_rows = _select_ready_rows(rows)
  if not ready_rows:
    raise SystemExit("[ERROR] no rows remain after reference-complex filtering")

  _write_ready_rows(args.ready_csv, ready_rows)
  print(
    f"[OK] ready index written: {args.ready_csv}\n"
    f"[INFO] rows_in={len(rows)} rows_ready={len(ready_rows)} complex_ok={ok_complexes} complex_fail={fail_complexes}"
  )

  extract_args = [
    "--dataset-csv", str(args.ready_csv),
    "--out-dir", str(args.epitope_dir),
    "--cutoff", str(args.cutoff),
  ]
  if not args.no_repair_target_pdbs:
    extract_args.append("--repair-target-pdbs")
  if not args.no_crop_targets:
    extract_args.extend(
      [
        "--crop-targets",
        "--cropped-target-dir", str(args.cropped_target_dir),
        "--crop-radius", str(args.crop_radius),
        "--crop-padding-residues", str(args.crop_padding_residues),
        "--crop-min-residues", str(args.crop_min_residues),
      ]
    )
  _run_python_step(data_prep_dir / "extract_epitopes_from_complexes.py", extract_args)

  repaired_ready_rows = _keep_rows_with_target_paths(_load_rows(args.ready_csv))
  if not repaired_ready_rows:
    raise SystemExit("[ERROR] no rows remain after target-structure repair")
  if len(repaired_ready_rows) != len(ready_rows):
    _write_ready_rows(args.ready_csv, repaired_ready_rows)
    print(f"[INFO] dropped_missing_targets={len(ready_rows) - len(repaired_ready_rows)}")

  _run_python_step(
    data_prep_dir / "build_model_inputs_native.py",
    [
      "--dataset-csv", str(args.ready_csv),
      "--out-root", str(args.native_root),
      "--models-root", str(args.models_root),
      "--epitope-dir", str(args.epitope_dir),
    ],
  )

  print("[OK] native input preparation complete")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())