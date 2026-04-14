#!/usr/bin/env python3
"""Select top-K candidates per (model, sample) from candidate_manifest.csv files.

For each model output directory, reads every sample's candidate_manifest.csv,
sorts by the model's internal score (or candidate_rank as fallback), takes the
top-K, and writes a new candidate_manifest.csv into an output directory tree
with the same structure.

Usage:
  python scripts/select_topk_manifest.py \\
    --input-root outputs/pilot_boltzgen \\
    --input-root outputs/pilot_rfantibody \\
    --output-root outputs/pilot_top10 \\
    --topk 10

Ranking strategy per model:
  boltzgen    : sort by native_score_value DESC (quality_score, higher is better)
  RFantibody  : no native score → sort by candidate_rank ASC (design index order)
  mber-open   : sort by native_score_value DESC (i_ptm, higher is better)
  germinal    : sort by native_score_value DESC (external_iptm / i_ptm, higher is better)
  default     : sort by candidate_rank ASC
"""
from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


# Models where higher native_score_value is better
SCORE_HIGHER_BETTER = {"boltzgen", "mber-open", "germinal"}
# Models where lower native_score_value is better
SCORE_LOWER_BETTER: set[str] = set()


def _sort_key_for_model(row: dict[str, str], model: str) -> tuple:
    """Return a sort key so that sorting ascending puts best candidates first."""
    raw_score = row.get("native_score_value", "").strip()
    try:
        score = float(raw_score)
    except ValueError:
        score = None

    rank_raw = row.get("candidate_rank", "").strip()
    try:
        rank = int(rank_raw)
    except ValueError:
        rank = 999999

    if score is not None:
        if model in SCORE_HIGHER_BETTER:
            return (0, -score, rank)   # negate so ascending = highest first
        if model in SCORE_LOWER_BETTER:
            return (0, score, rank)    # ascending = lowest first
        # Unknown model with a score: treat higher as better
        return (0, -score, rank)
    else:
        # No score: fall back to candidate_rank
        return (1, rank, 0)


def select_topk(
    input_roots: list[Path],
    output_root: Path,
    topk: int = 10,
) -> dict[str, int]:
    """Process all manifests and write filtered ones.

    Returns summary dict {model/sample_id: n_selected}.
    """
    summary: dict[str, int] = {}

    for root in input_roots:
        if not root.is_dir():
            print(f"[WARN] input-root does not exist: {root}")
            continue

        # Expected layout: root/<model>/<sample_id>/candidate_manifest.csv
        for model_dir in sorted(root.iterdir()):
            if not model_dir.is_dir():
                continue
            model = model_dir.name
            if model == "manifest.csv":
                continue

            for sample_dir in sorted(model_dir.iterdir()):
                if not sample_dir.is_dir():
                    continue
                manifest_path = sample_dir / "candidate_manifest.csv"
                if not manifest_path.exists():
                    print(f"[SKIP] no manifest: {manifest_path}")
                    continue

                sample_id = sample_dir.name

                # Read all ok rows
                with manifest_path.open("r", encoding="utf-8", newline="") as f:
                    reader = csv.DictReader(f)
                    fieldnames = list(reader.fieldnames or [])
                    all_rows = list(reader)

                ok_rows = [r for r in all_rows if r.get("status") == "ok"]
                failed_rows = [r for r in all_rows if r.get("status") != "ok"]

                if not ok_rows:
                    print(f"[WARN] no ok candidates: {model}/{sample_id}")
                    # Still write empty manifest so downstream knows this was attempted
                    top_rows: list[dict[str, str]] = failed_rows[:1]
                else:
                    # Sort and select top-K
                    sorted_ok = sorted(ok_rows, key=lambda r: _sort_key_for_model(r, model))
                    top_rows = sorted_ok[:topk]

                # Write to output
                out_dir = output_root / model / sample_id
                out_dir.mkdir(parents=True, exist_ok=True)
                out_manifest = out_dir / "candidate_manifest.csv"

                with out_manifest.open("w", encoding="utf-8", newline="") as f:
                    w = csv.DictWriter(f, fieldnames=fieldnames)
                    w.writeheader()
                    w.writerows(top_rows)

                key = f"{model}/{sample_id}"
                summary[key] = len(top_rows)
                print(
                    f"[OK] {key}: {len(ok_rows)} ok → top {len(top_rows)} written"
                    f"  (score={top_rows[0].get('native_score_value','N/A')} "
                    f"rank={top_rows[0].get('candidate_rank','?')})"
                    if top_rows else f"[OK] {key}: 0 candidates"
                )

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Select top-K candidates from model output manifests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input-root", action="append", type=Path, required=True,
        dest="input_roots",
        help="Model output directory (<model>/<sample_id>/candidate_manifest.csv). Repeatable.",
    )
    parser.add_argument(
        "--output-root", type=Path, required=True,
        help="Destination directory. Same layout as input-root but with top-K manifests.",
    )
    parser.add_argument(
        "--topk", type=int, default=10,
        help="Number of candidates to keep per (model, sample). Default: 10.",
    )
    args = parser.parse_args()

    print(f"[select_topk] top-K={args.topk}")
    print(f"[select_topk] input_roots={[str(r) for r in args.input_roots]}")
    print(f"[select_topk] output_root={args.output_root}")

    summary = select_topk(
        input_roots=args.input_roots,
        output_root=args.output_root,
        topk=args.topk,
    )

    print(f"\n[Done] {len(summary)} (model, sample) pairs processed.")
    total = sum(summary.values())
    print(f"[Done] Total candidates selected: {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
