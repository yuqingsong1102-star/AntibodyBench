#!/usr/bin/env python3
"""Aggregate: merge results from the three evaluation steps into a single CSV.

Input:  chai1_results.csv, interface_results.csv, sequence_results.csv
Output: eval_merged.csv — one row per candidate with all metrics joined.

All metrics are continuous. No pass/fail gates.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


# The merged output keeps all fields from all steps, with prefixed namespacing.
# Common fields (candidate_id, sample_id, model) are shared.
COMMON_FIELDS = ["candidate_id", "sample_id", "model"]


def run_aggregate(
    *,
    chai1_csv: Path,
    interface_csv: Path,
    sequence_csv: Path,
    out_csv: Path,
) -> Path:
    """Merge three step CSVs into one unified result table."""
    print("[Aggregate] Merging evaluation results...")

    # Load each step
    chai1_rows = {r["candidate_id"]: r for r in _read_csv(chai1_csv) if r.get("candidate_id")}
    iface_rows = {r["candidate_id"]: r for r in _read_csv(interface_csv) if r.get("candidate_id")}
    seq_rows = {r["candidate_id"]: r for r in _read_csv(sequence_csv) if r.get("candidate_id")}

    # Collect all candidate IDs (union)
    all_ids = sorted(set(chai1_rows) | set(iface_rows) | set(seq_rows))

    if not all_ids:
        print("[Aggregate] No candidates found in any step CSV.")
        return out_csv

    # Determine output field order
    # Common first, then chai1, interface, sequence (excluding duplicates)
    chai1_fields = []
    iface_fields = []
    seq_fields = []

    if chai1_rows:
        sample = next(iter(chai1_rows.values()))
        chai1_fields = [k for k in sample if k not in COMMON_FIELDS]
    if iface_rows:
        sample = next(iter(iface_rows.values()))
        iface_fields = [k for k in sample if k not in COMMON_FIELDS and k not in chai1_fields]
    if seq_rows:
        sample = next(iter(seq_rows.values()))
        # Skip pred_sequence if already in chai1 (avoid duplication)
        existing = set(COMMON_FIELDS) | set(chai1_fields) | set(iface_fields)
        seq_fields = [k for k in sample if k not in existing]

    all_fields = COMMON_FIELDS + chai1_fields + iface_fields + seq_fields

    # Merge
    merged_rows: list[dict[str, str]] = []
    for cid in all_ids:
        row: dict[str, str] = {k: "" for k in all_fields}
        row["candidate_id"] = cid

        # Fill from each source (later sources don't overwrite non-empty values)
        for source in [chai1_rows.get(cid, {}), iface_rows.get(cid, {}), seq_rows.get(cid, {})]:
            for k, v in source.items():
                if k in row and (not row[k] or row[k] == ""):
                    row[k] = v

        merged_rows.append(row)

    # Write
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=all_fields)
        w.writeheader()
        w.writerows(merged_rows)

    # Summary
    n_ok_chai = sum(1 for r in merged_rows if r.get("chai1_status") == "ok")
    n_ok_iface = sum(1 for r in merged_rows if r.get("interface_status") == "ok")
    n_ok_seq = sum(1 for r in merged_rows if r.get("sequence_status") == "ok")

    print(f"[Aggregate] Merged {len(merged_rows)} candidates")
    print(f"[Aggregate]   Chai-1 OK: {n_ok_chai}")
    print(f"[Aggregate]   Interface OK: {n_ok_iface}")
    print(f"[Aggregate]   Sequence OK: {n_ok_seq}")
    print(f"[Aggregate] Written to {out_csv}")
    return out_csv


def merge_extra_csv(
    *,
    base_csv: Path,
    extra_csv: Path,
    out_csv: Path,
) -> Path:
    """Merge an extra CSV (e.g. diversity_epitope) into a base CSV by candidate_id."""
    if not base_csv.exists():
        print(f"[Aggregate] Base CSV not found: {base_csv}")
        return out_csv

    base_rows = _read_csv(base_csv)
    extra_map: dict[str, dict[str, str]] = {}
    extra_fields: list[str] = []
    if extra_csv.exists():
        extra_rows = _read_csv(extra_csv)
        if extra_rows:
            extra_fields = [k for k in extra_rows[0] if k not in COMMON_FIELDS]
            extra_map = {r["candidate_id"]: r for r in extra_rows if r.get("candidate_id")}

    if not extra_map:
        # Nothing to merge, just copy base
        if base_csv != out_csv:
            import shutil
            shutil.copy2(base_csv, out_csv)
        return out_csv

    # Build merged fieldset
    base_fields = list(base_rows[0].keys()) if base_rows else COMMON_FIELDS
    new_fields = [f for f in extra_fields if f not in base_fields]
    all_fields = base_fields + new_fields

    merged = []
    for r in base_rows:
        cid = r.get("candidate_id", "")
        row = {k: r.get(k, "") for k in all_fields}
        if cid in extra_map:
            for k in new_fields:
                row[k] = extra_map[cid].get(k, "")
        merged.append(row)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=all_fields)
        w.writeheader()
        w.writerows(merged)

    print(f"[Aggregate] Merged {len(new_fields)} extra fields into {out_csv}")
    return out_csv


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate evaluation step results")
    parser.add_argument(
        "--chai1-csv", type=Path,
        default=Path("outputs/evaluation/results/chai1_results.csv"),
    )
    parser.add_argument(
        "--interface-csv", type=Path,
        default=Path("outputs/evaluation/results/interface_results.csv"),
    )
    parser.add_argument(
        "--sequence-csv", type=Path,
        default=Path("outputs/evaluation/results/sequence_results.csv"),
    )
    parser.add_argument(
        "--out-csv", type=Path,
        default=Path("outputs/evaluation/results/eval_merged.csv"),
    )
    args = parser.parse_args()

    run_aggregate(
        chai1_csv=args.chai1_csv,
        interface_csv=args.interface_csv,
        sequence_csv=args.sequence_csv,
        out_csv=args.out_csv,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
