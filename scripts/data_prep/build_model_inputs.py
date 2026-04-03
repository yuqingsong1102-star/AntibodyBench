#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


SUPPORTED_MODELS = ("RFantibody", "germinal", "boltzgen", "BindCraft")


def _run_preprocess(
  *,
  root: Path,
  model: str,
  dataset_csv: Path,
  af3_input: Path,
  out_dir: Path,
) -> None:
  preprocess = root / "algorithms" / model / "preprocess.py"
  if not preprocess.exists():
    raise FileNotFoundError(f"preprocess not found: {preprocess}")

  cmd = [
    sys.executable,
    str(preprocess),
    "--af3-input",
    str(af3_input),
    "--dataset-csv",
    str(dataset_csv),
    "--out-dir",
    str(out_dir),
  ]
  print(f"[RUN] {' '.join(cmd)}")
  subprocess.run(cmd, check=True)


def main() -> int:
  parser = argparse.ArgumentParser(description="Build per-model input files under data/model_inputs/")
  parser.add_argument(
    "--dataset-csv",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "inputs" / "dataset_index_ready.csv",
    help="Dataset index used to generate model inputs.",
  )
  parser.add_argument(
    "--af3-input",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "inputs" / "alphafold3_inputs.json",
    help="AF3-style optional metadata input.",
  )
  parser.add_argument(
    "--out-root",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "model_inputs",
    help="Output root directory for per-model inputs.",
  )
  parser.add_argument(
    "--models",
    nargs="+",
    choices=SUPPORTED_MODELS,
    default=list(SUPPORTED_MODELS),
    help="Models to prepare.",
  )
  parser.add_argument(
    "--clean",
    action="store_true",
    help="Delete existing per-model output dirs before regeneration.",
  )
  args = parser.parse_args()

  root = Path(__file__).resolve().parents[2]

  if not args.dataset_csv.exists():
    raise SystemExit(f"[ERROR] dataset csv not found: {args.dataset_csv}")
  if not args.af3_input.exists():
    raise SystemExit(f"[ERROR] af3 input not found: {args.af3_input}")

  args.out_root.mkdir(parents=True, exist_ok=True)

  for model in args.models:
    out_dir = args.out_root / model
    if args.clean and out_dir.exists():
      shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _run_preprocess(
      root=root,
      model=model,
      dataset_csv=args.dataset_csv,
      af3_input=args.af3_input,
      out_dir=out_dir,
    )
    print(f"[OK] model inputs ready: {out_dir}")

  print("[DONE] All requested model inputs generated.")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

