#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import urllib.error
import urllib.request
from pathlib import Path


def _download(url: str, out_path: Path, timeout_sec: int) -> bool:
  try:
    with urllib.request.urlopen(url, timeout=timeout_sec) as resp:
      data = resp.read()
    if not data:
      return False
    out_path.write_bytes(data)
    return True
  except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
    return False


def _is_nonempty_file(path: Path) -> bool:
  return path.exists() and path.is_file() and path.stat().st_size > 20


def _base_pdb_id(raw: str) -> str:
  return (raw or "").strip().split("-")[0].lower()


def main() -> int:
  parser = argparse.ArgumentParser(description="Download full complex structures and fill reference_complex_path")
  parser.add_argument(
    "--input-csv",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "inputs" / "dataset_index.csv",
    help="Input dataset index CSV.",
  )
  parser.add_argument(
    "--output-csv",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "inputs" / "dataset_index.csv",
    help="Output CSV (can be same as input).",
  )
  parser.add_argument(
    "--complex-dir",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "raw" / "reference_complexes",
    help="Directory to store downloaded full-complex structures.",
  )
  parser.add_argument(
    "--timeout-sec",
    type=int,
    default=8,
    help="HTTP timeout seconds per request (default: 8).",
  )
  args = parser.parse_args()

  with args.input_csv.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))
    if not rows:
      raise SystemExit("input csv is empty")
    fieldnames = list(rows[0].keys())

  for k in ("reference_complex_path", "reference_complex_status", "reference_complex_note"):
    if k not in fieldnames:
      fieldnames.append(k)

  args.complex_dir.mkdir(parents=True, exist_ok=True)

  ok = 0
  missing_pdb = 0
  download_fail = 0
  cache: dict[str, tuple[str, str]] = {}
  for row in rows:
    pdb_id = _base_pdb_id(str(row.get("pdb_id", "")))
    if not pdb_id:
      row["reference_complex_status"] = "fail"
      row["reference_complex_note"] = "missing_pdb_id"
      missing_pdb += 1
      continue

    if pdb_id in cache:
      status, value = cache[pdb_id]
      if status == "ok":
        row["reference_complex_path"] = value
        row["reference_complex_status"] = "ok"
        row["reference_complex_note"] = ""
        ok += 1
      else:
        row["reference_complex_status"] = "fail"
        row["reference_complex_note"] = value
      continue

    out_pdb = args.complex_dir / f"{pdb_id}.pdb"
    out_cif = args.complex_dir / f"{pdb_id}.cif"
    if _is_nonempty_file(out_pdb):
      cache[pdb_id] = ("ok", str(out_pdb.resolve()))
      row["reference_complex_path"] = str(out_pdb.resolve())
      row["reference_complex_status"] = "ok"
      row["reference_complex_note"] = ""
      ok += 1
      continue
    if _is_nonempty_file(out_cif):
      cache[pdb_id] = ("ok", str(out_cif.resolve()))
      row["reference_complex_path"] = str(out_cif.resolve())
      row["reference_complex_status"] = "ok"
      row["reference_complex_note"] = ""
      ok += 1
      continue

    pdb_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    if _download(pdb_url, out_pdb, timeout_sec=args.timeout_sec) and _is_nonempty_file(out_pdb):
      cache[pdb_id] = ("ok", str(out_pdb.resolve()))
      row["reference_complex_path"] = str(out_pdb.resolve())
      row["reference_complex_status"] = "ok"
      row["reference_complex_note"] = ""
      ok += 1
      continue

    cif_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    if _download(cif_url, out_cif, timeout_sec=args.timeout_sec) and _is_nonempty_file(out_cif):
      cache[pdb_id] = ("ok", str(out_cif.resolve()))
      row["reference_complex_path"] = str(out_cif.resolve())
      row["reference_complex_status"] = "ok"
      row["reference_complex_note"] = "downloaded_cif"
      ok += 1
      continue

    cache[pdb_id] = ("fail", "download_failed")
    row["reference_complex_status"] = "fail"
    row["reference_complex_note"] = "download_failed"
    download_fail += 1

  args.output_csv.parent.mkdir(parents=True, exist_ok=True)
  with args.output_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

  print(
    f"[OK] wrote: {args.output_csv}\n"
    f"[INFO] rows={len(rows)} ok={ok} missing_pdb={missing_pdb} download_fail={download_fail}\n"
    f"[INFO] complex_dir={args.complex_dir}"
  )
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

