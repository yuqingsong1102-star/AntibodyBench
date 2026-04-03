#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


MODELS = ("RFantibody", "germinal", "boltzgen", "BindCraft")


def _safe(s: str) -> str:
  return (s or "").strip()


def _write_text(path: Path, text: str) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_text(text, encoding="utf-8")


def _build_rfantibody(sample_dir: Path, row: dict[str, str], models_root: Path) -> list[Path]:
  sample_id = _safe(row.get("sample_id", "sample"))
  antigen_chain = _safe(row.get("antigen_chain", "A"))
  target_pdb = _safe(row.get("reference_structure_path", ""))
  framework_default = models_root / "RFantibody" / "scripts" / "examples" / "example_inputs" / "h-NbBCII10.pdb"

  files: list[Path] = []
  cfg = sample_dir / "rfantibody_config.env"
  _write_text(
    cfg,
    "\n".join(
      [
        f'SAMPLE_ID="{sample_id}"',
        f'TARGET_PDB="{target_pdb}"',
        f'FRAMEWORK_PDB="{framework_default}"',
        f'HOTSPOTS="{antigen_chain}__FILL_ME__"',
        'DESIGN_LOOPS="H1:7,H2:6,H3:5-13"',
        "NUM_DESIGNS=10",
        "NUM_SEQS=4",
        "NUM_RECYCLES=10",
        "",
      ]
    ),
  )
  files.append(cfg)
  return files


def _build_germinal(sample_dir: Path, row: dict[str, str]) -> list[Path]:
  target_pdb = _safe(row.get("reference_structure_path", ""))
  antigen_chain = _safe(row.get("antigen_chain", "A"))
  # Germinal convention: binder_chain is usually B in target config.
  # Here we preserve dataset value and leave hotspots placeholder.
  binder_chain = _safe(row.get("antibody_chain", "B"))
  sample_id = _safe(row.get("sample_id", "sample"))

  files: list[Path] = []
  yaml_path = sample_dir / "target.yaml"
  _write_text(
    yaml_path,
    "\n".join(
      [
        f'target_name: "{sample_id}"',
        f'target_pdb_path: "{target_pdb}"',
        f'target_chain: "{antigen_chain}"',
        f'binder_chain: "{binder_chain}"',
        f'target_hotspots: "{antigen_chain}__FILL_ME__"',
        "",
      ]
    ),
  )
  files.append(yaml_path)

  overrides = sample_dir / "run_overrides.txt"
  _write_text(
    overrides,
    "\n".join(
      [
        f"target.target_pdb_path='{target_pdb}'",
        f"target.target_chain='{antigen_chain}'",
        f"target.binder_chain='{binder_chain}'",
        f"target.target_hotspots='{antigen_chain}__FILL_ME__'",
        "",
      ]
    ),
  )
  files.append(overrides)
  return files


def _build_bindcraft(sample_dir: Path, row: dict[str, str]) -> list[Path]:
  sample_id = _safe(row.get("sample_id", "sample"))
  target_pdb = _safe(row.get("reference_structure_path", ""))
  antigen_chain = _safe(row.get("antigen_chain", "A"))

  files: list[Path] = []
  settings = sample_dir / "settings_target.json"
  payload = {
    "design_path": str((sample_dir / "designs").resolve()),
    "binder_name": sample_id,
    "starting_pdb": target_pdb,
    "chains": antigen_chain,
    "target_hotspot_residues": None,  # Fill with e.g. "A37,A39,A41"
    "lengths": [20, 35],
    "number_of_final_designs": 10,
  }
  _write_text(settings, json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
  files.append(settings)
  return files


def _build_boltzgen(sample_dir: Path, row: dict[str, str], models_root: Path) -> list[Path]:
  target_path = _safe(row.get("reference_structure_path", ""))
  antigen_chain = _safe(row.get("antigen_chain", "A"))
  # Use bundled scaffold examples from boltzgen repo as a starting point.
  scaffold_candidates = [
    models_root / "boltzgen" / "example" / "nanobody_scaffolds" / "7eow.yaml",
    models_root / "boltzgen" / "example" / "nanobody_scaffolds" / "7xl0.yaml",
  ]

  files: list[Path] = []
  yaml_path = sample_dir / "design_spec.yaml"
  lines = [
    "entities:",
    "  - file:",
    f'      path: "{target_path}"',
    "      include:",
    "        - chain:",
    f'            id: "{antigen_chain}"',
    "  - file:",
    "      path:",
  ]
  for cand in scaffold_candidates:
    lines.append(f'        - "{cand}"')
  lines += [
    "",
    "# Optional: add hotspot-focused constraints / include-exclude masks per boltzgen examples.",
    "",
  ]
  _write_text(yaml_path, "\n".join(lines))
  files.append(yaml_path)
  return files


def main() -> int:
  parser = argparse.ArgumentParser(description="Generate model-native input templates per sample.")
  parser.add_argument(
    "--dataset-csv",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "inputs" / "dataset_index_ready.csv",
    help="Input benchmark index CSV.",
  )
  parser.add_argument(
    "--out-root",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "model_inputs_native",
    help="Output root for native input templates.",
  )
  parser.add_argument(
    "--models-root",
    type=Path,
    default=Path(__file__).resolve().parents[3] / "models",
    help="Root path of original model repositories.",
  )
  args = parser.parse_args()

  if not args.dataset_csv.exists():
    raise SystemExit(f"[ERROR] dataset csv not found: {args.dataset_csv}")

  with args.dataset_csv.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))

  args.out_root.mkdir(parents=True, exist_ok=True)
  manifest_rows: list[dict[str, str]] = []

  for row in rows:
    sample_id = _safe(row.get("sample_id", ""))
    if not sample_id:
      continue

    for model in MODELS:
      sample_dir = args.out_root / model / sample_id
      sample_dir.mkdir(parents=True, exist_ok=True)

      if model == "RFantibody":
        files = _build_rfantibody(sample_dir, row, args.models_root)
      elif model == "germinal":
        files = _build_germinal(sample_dir, row)
      elif model == "BindCraft":
        files = _build_bindcraft(sample_dir, row)
      elif model == "boltzgen":
        files = _build_boltzgen(sample_dir, row, args.models_root)
      else:
        files = []

      for fp in files:
        manifest_rows.append({"model": model, "sample_id": sample_id, "path": str(fp.resolve())})

  manifest = args.out_root / "_native_manifest.csv"
  with manifest.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["model", "sample_id", "path"])
    writer.writeheader()
    writer.writerows(manifest_rows)

  print(f"[OK] native model inputs generated: {args.out_root}")
  print(f"[INFO] manifest: {manifest}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

