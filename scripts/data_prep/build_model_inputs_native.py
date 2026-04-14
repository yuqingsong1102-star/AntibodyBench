#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


MODELS = ("RFantibody", "germinal", "boltzgen", "mber-open")

DEFAULT_MBER_MASKED_VHH = (
  "EVQLVESGGGLVQPGGSLRLSCAASG*********WFRQAPGKEREF***********"
  "NADSVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYC************WGQGTLVTVSS"
)


def _safe(s: str) -> str:
  return (s or "").strip()


def _write_text(path: Path, text: str) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_text(text, encoding="utf-8")


def _read_json(path: Path) -> dict:
  if not path.exists():
    return {}
  try:
    return json.loads(path.read_text(encoding="utf-8"))
  except Exception:
    return {}


def _has_af_params(path: Path) -> bool:
  if not path.exists() or not path.is_dir():
    return False
  return any(path.glob("params_model_*.npz"))


def _resolve_germinal_af_params_dir(models_root: Path) -> Path:
  candidates = [models_root / "germinal" / "params"]
  for candidate in candidates:
    if _has_af_params(candidate):
      return candidate
  return candidates[0]


def _get_hotspots(epitope: dict, antigen_chain: str, *, allow_placeholder: bool) -> str | None:
  hotspot_string = _safe(str(epitope.get("hotspot_string") or ""))
  hotspot_count = int(epitope.get("hotspot_count") or 0)
  if _safe(str(epitope.get("status") or "")) == "ok" and hotspot_count > 0 and hotspot_string:
    return hotspot_string
  if not allow_placeholder:
    return None
  return f"{antigen_chain}__FILL_ME__"


def _split_chains(chain_field: str) -> list[str]:
  raw = _safe(chain_field)
  if not raw:
    return []
  for sep in (";", ",", " "):
    if sep in raw:
      return [part.strip() for part in raw.split(sep) if part.strip()]
  return [raw]


def _pick_designed_binder_chain(row: dict[str, str]) -> str:
  antigen_chains = set(_split_chains(_safe(row.get("antigen_chain", ""))))
  preferred = ["B", "H", "L", "C", "D", "E", "F", "G", "I", "J", "K", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "A"]
  for chain_id in preferred:
    if chain_id not in antigen_chains:
      return chain_id
  return "B"


def _get_target_pdb(row: dict[str, str], epitope: dict) -> str:
  return (
    _safe(str(row.get("design_structure_path") or ""))
    or _safe(str(epitope.get("design_structure_path") or ""))
    or _safe(str(row.get("reference_structure_path") or ""))
  )


def _build_rfantibody(sample_dir: Path, row: dict[str, str], models_root: Path, epitope: dict) -> list[Path]:
  sample_id = _safe(row.get("sample_id", "sample"))
  antigen_chain = _safe(row.get("antigen_chain", "A"))
  target_pdb = _get_target_pdb(row, epitope)
  framework_default = models_root / "RFantibody" / "scripts" / "examples" / "example_inputs" / "h-NbBCII10.pdb"
  hotspots = _get_hotspots(epitope, antigen_chain, allow_placeholder=True) or f"{antigen_chain}__FILL_ME__"

  files: list[Path] = []
  cfg = sample_dir / "rfantibody_config.env"
  _write_text(
    cfg,
    "\n".join(
      [
        f'SAMPLE_ID="{sample_id}"',
        f'TARGET_PDB="{target_pdb}"',
        f'FRAMEWORK_PDB="{framework_default}"',
        f'HOTSPOTS="{hotspots}"',
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


def _build_germinal(sample_dir: Path, row: dict[str, str], models_root: Path, epitope: dict) -> list[Path]:
  target_pdb = _get_target_pdb(row, epitope)
  antigen_chain = _safe(row.get("antigen_chain", "A"))
  binder_chain = _pick_designed_binder_chain(row)
  hotspots = _get_hotspots(epitope, antigen_chain, allow_placeholder=True) or f"{antigen_chain}__FILL_ME__"
  af_params_dir = _resolve_germinal_af_params_dir(models_root)
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
        f'target_hotspots: "{hotspots}"',
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
        f"target.target_name='{sample_id}'",
        f"target.target_pdb_path='{target_pdb}'",
        f"target.target_chain='{antigen_chain}'",
        f"target.binder_chain='{binder_chain}'",
        f"target.target_hotspots='{hotspots}'",
        "target.hotspot_residue=null",
        f"af_params_dir='{af_params_dir}'",
        "",
      ]
    ),
  )
  files.append(overrides)
  return files


def _build_mber_open(sample_dir: Path, row: dict[str, str], epitope: dict) -> list[Path]:
  sample_id = _safe(row.get("sample_id", "sample"))
  target_pdb = _get_target_pdb(row, epitope)
  antigen_chain = _safe(row.get("antigen_chain", "A"))
  hotspots = _get_hotspots(epitope, antigen_chain, allow_placeholder=False)

  files: list[Path] = []
  settings = sample_dir / "mber_vhh_settings.json"
  payload = {
    "output": {
      "dir": "__OUTPUT_DIR__",
      "skip_animations": True,
      "skip_pickle": True,
      "skip_png": True,
    },
    "target": {
      "pdb": target_pdb,
      "name": sample_id,
      "chains": antigen_chain,
      "hotspots": [part for part in (hotspots or "").split(",") if part],
    },
    "binder": {
      "masked_sequence": DEFAULT_MBER_MASKED_VHH,
    },
    "stopping": {
      "num_accepted": 3,
      "max_trajectories": 50,
    },
    "filters": {
      "min_iptm": 0.75,
      "min_plddt": 0.70,
    },
  }
  _write_text(settings, json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
  files.append(settings)
  return files


def _build_boltzgen(sample_dir: Path, row: dict[str, str], models_root: Path, epitope: dict) -> list[Path]:
  target_path = _get_target_pdb(row, epitope)
  antigen_chain = _safe(row.get("antigen_chain", "A"))
  # Use bundled scaffold examples from boltzgen repo as a starting point.
  scaffold_candidates = [
    models_root / "boltzgen" / "example" / "nanobody_scaffolds" / "7eow.yaml",
    models_root / "boltzgen" / "example" / "nanobody_scaffolds" / "7xl0.yaml",
  ]
  hotspots = _get_hotspots(epitope, antigen_chain, allow_placeholder=False)

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
    f'# target_hotspots_from_epitopes: "{hotspots or ""}"',
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
    default=Path(__file__).resolve().parents[2] / "data" / "dataset_index_ready.csv",
    help="Input benchmark index CSV.",
  )
  parser.add_argument(
    "--out-root",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "native_inputs",
    help="Output root for native input templates.",
  )
  parser.add_argument(
    "--models-root",
    type=Path,
    default=Path(__file__).resolve().parents[3] / "models",
    help="Root path of original model repositories.",
  )
  parser.add_argument(
    "--epitope-dir",
    type=Path,
    default=Path(__file__).resolve().parents[2] / "data" / "epitopes",
    help="Directory of per-sample epitope json files.",
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
    epitope = _read_json(args.epitope_dir / f"{sample_id}.json")

    for model in MODELS:
      sample_dir = args.out_root / model / sample_id
      sample_dir.mkdir(parents=True, exist_ok=True)

      if model == "RFantibody":
        files = _build_rfantibody(sample_dir, row, args.models_root, epitope)
      elif model == "germinal":
        files = _build_germinal(sample_dir, row, args.models_root, epitope)
      elif model == "mber-open":
        files = _build_mber_open(sample_dir, row, epitope)
      elif model == "boltzgen":
        files = _build_boltzgen(sample_dir, row, args.models_root, epitope)
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

