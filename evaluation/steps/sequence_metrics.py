#!/usr/bin/env python3
"""Step 3: Sequence-level analysis.

Pure sequence analysis — no GPU, no structure prediction.
Computes developability flags, physicochemical properties, CDR annotation,
and sequence-level quality indicators.

Input:  chai1_results.csv (for sequences) or model candidate manifests
Output: sequence_results.csv with continuous sequence metrics per candidate.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

from evaluation.structure_inference import parse_structure, protein_chain_sequences

# ---------------------------------------------------------------------------
# Output schema
# ---------------------------------------------------------------------------
SEQUENCE_FIELDS = [
    "candidate_id",
    "sample_id",
    "model",
    "sequence_status",
    "pred_sequence",
    "seq_length",
    # Physicochemical
    "net_charge",
    "isoelectric_point_approx",
    "molecular_weight_approx",
    "gravy",  # grand average of hydropathicity
    "sequence_entropy",
    # Liability motifs
    "free_cys_count",
    "glycosylation_motif_count",
    "deamidation_motif_count",
    "isomerization_motif_count",
    "total_liability_count",
    # Composition
    "aromatic_fraction",
    "charged_fraction",
    "hydrophobic_fraction",
    # VHH-specific
    "has_vhh_hallmarks",
    "framework_gapped_positions",
    # CDR lengths (Chothia-like heuristic)
    "cdr1_length",
    "cdr2_length",
    "cdr3_length",
    "cdr3_sequence",
]

# Physicochemical constants
AA_MW = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.16,
    "E": 147.13, "Q": 146.15, "G": 75.03, "H": 155.16, "I": 131.17,
    "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
    "S": 105.09, "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
}
AA_CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1}
AA_HYDROPATHY = {  # Kyte-Doolittle
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "E": -3.5, "Q": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
AROMATIC = set("FWY")
CHARGED = set("DEKRH")
HYDROPHOBIC = set("AVILMFWP")

GLYCO_RE = re.compile(r"N[^P][ST]")
DEAMIDATION_RE = re.compile(r"N[GS]")
ISOMERIZATION_RE = re.compile(r"D[GST]")

# VHH hallmarks (positions in Kabat numbering, approximate in sequence)
# Position 42: F/Y (framework 2), Position 49: E/Q, Position 50: R/C
# Simplified: check for characteristic VHH substitutions in FR2 region
VHH_FR2_PATTERN = re.compile(
    r"[FY][A-Z]{6,8}[EQ][A-Z][RC]",
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _write_results(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SEQUENCE_FIELDS)
        w.writeheader()
        w.writerows(rows)


def _r(value: float | None, digits: int = 4) -> str:
    return "" if value is None else str(round(value, digits))


def _motif_count(pattern: re.Pattern, seq: str) -> int:
    return sum(1 for _ in pattern.finditer(seq))


def _sequence_entropy(seq: str) -> float | None:
    if not seq:
        return None
    counts: dict[str, int] = {}
    for aa in seq:
        counts[aa] = counts.get(aa, 0) + 1
    total = len(seq)
    return -sum((c / total) * math.log(c / total, 2) for c in counts.values())


def _net_charge(seq: str) -> float:
    return sum(AA_CHARGE.get(aa, 0.0) for aa in seq)


def _approx_pi(seq: str) -> float | None:
    """Very rough isoelectric point estimate."""
    if not seq:
        return None
    # Count ionizable residues
    n_asp = seq.count("D")
    n_glu = seq.count("E")
    n_lys = seq.count("K")
    n_arg = seq.count("R")
    n_his = seq.count("H")
    # Simplified: pI ≈ (sum of pKa of basic) / n_basic if net positive at pH 7
    # This is a very rough approximation
    positive = n_lys + n_arg + n_his
    negative = n_asp + n_glu
    if positive + negative == 0:
        return 7.0
    # Rough: pI shifts higher with more basic, lower with more acidic
    return 7.0 + 0.5 * (positive - negative) / max(1, positive + negative) * 4.0


def _molecular_weight(seq: str) -> float:
    """Approximate molecular weight in Daltons."""
    return sum(AA_MW.get(aa, 120.0) for aa in seq) - (len(seq) - 1) * 18.015


def _gravy(seq: str) -> float | None:
    if not seq:
        return None
    return sum(AA_HYDROPATHY.get(aa, 0.0) for aa in seq) / len(seq)


def _estimate_cdr_lengths(seq: str) -> tuple[int, int, int, str]:
    """Heuristic CDR annotation for VHH/VH sequences.

    Uses a simplified positional heuristic based on typical VHH framework lengths:
    - FR1: ~25 aa
    - CDR1: ~6-12 aa (positions 26-35ish)
    - FR2: ~17 aa
    - CDR2: ~3-9 aa (positions 52-58ish)
    - FR3: ~32 aa
    - CDR3: variable (positions 95-102ish)
    - FR4: ~11 aa

    Returns (cdr1_len, cdr2_len, cdr3_len, cdr3_seq).
    """
    n = len(seq)
    if n < 80:
        return 0, 0, 0, ""

    # Approximate positions (0-indexed)
    cdr1_start = 25
    cdr1_end = min(35, n - 60)
    if cdr1_end <= cdr1_start:
        cdr1_end = cdr1_start + 6

    cdr2_start = cdr1_end + 15  # FR2 ~ 15-17 aa
    cdr2_end = cdr2_start + 7
    if cdr2_end > n - 30:
        cdr2_end = cdr2_start + 3

    # CDR3: typically starts around position 95-100 and ends ~11 from C-terminus
    fr4_len = 11
    cdr3_end = max(cdr2_end + 20, n - fr4_len)
    cdr3_start = max(cdr2_end + 30, cdr3_end - 25)
    if cdr3_start >= cdr3_end:
        cdr3_start = cdr3_end - 5

    cdr1_len = max(0, cdr1_end - cdr1_start)
    cdr2_len = max(0, cdr2_end - cdr2_start)
    cdr3_len = max(0, cdr3_end - cdr3_start)
    cdr3_seq = seq[cdr3_start:cdr3_end] if cdr3_start < n else ""

    return cdr1_len, cdr2_len, cdr3_len, cdr3_seq


def _analyze_sequence(seq: str) -> dict[str, str]:
    """Analyze a single VHH candidate sequence."""
    if not seq:
        return {k: "" for k in SEQUENCE_FIELDS if k not in ("candidate_id", "sample_id", "model")}

    n = len(seq)
    cys_count = seq.count("C")
    glyco = _motif_count(GLYCO_RE, seq)
    deamid = _motif_count(DEAMIDATION_RE, seq)
    isom = _motif_count(ISOMERIZATION_RE, seq)
    total_liability = glyco + deamid + isom

    aromatic_frac = sum(1 for aa in seq if aa in AROMATIC) / n
    charged_frac = sum(1 for aa in seq if aa in CHARGED) / n
    hydro_frac = sum(1 for aa in seq if aa in HYDROPHOBIC) / n

    # VHH hallmarks
    has_vhh = bool(VHH_FR2_PATTERN.search(seq[30:65])) if n > 65 else False

    cdr1, cdr2, cdr3, cdr3_seq = _estimate_cdr_lengths(seq)

    return {
        "pred_sequence": seq,
        "seq_length": str(n),
        "net_charge": _r(_net_charge(seq)),
        "isoelectric_point_approx": _r(_approx_pi(seq)),
        "molecular_weight_approx": _r(_molecular_weight(seq), 1),
        "gravy": _r(_gravy(seq)),
        "sequence_entropy": _r(_sequence_entropy(seq)),
        "free_cys_count": str(cys_count),
        "glycosylation_motif_count": str(glyco),
        "deamidation_motif_count": str(deamid),
        "isomerization_motif_count": str(isom),
        "total_liability_count": str(total_liability),
        "aromatic_fraction": _r(aromatic_frac),
        "charged_fraction": _r(charged_frac),
        "hydrophobic_fraction": _r(hydro_frac),
        "has_vhh_hallmarks": "1" if has_vhh else "0",
        "framework_gapped_positions": "",
        "cdr1_length": str(cdr1),
        "cdr2_length": str(cdr2),
        "cdr3_length": str(cdr3),
        "cdr3_sequence": cdr3_seq,
    }


def _extract_seq_from_structure(pdb_path: Path, chain_hint: str = "") -> str:
    """Extract sequence from PDB when no explicit sequence is available."""
    try:
        structure = parse_structure(pdb_path)
        seqs = protein_chain_sequences(structure)
        if chain_hint and chain_hint in seqs:
            return seqs[chain_hint]
        if seqs:
            return min(seqs.values(), key=len)
    except Exception:
        pass
    return ""


def run_step3_sequence(
    *,
    chai1_csv: Path | None = None,
    native_roots: list[Path] | None = None,
    index_csv: Path,
    out_csv: Path,
    models: list[str] | None = None,
) -> Path:
    """Run sequence analysis on all candidates.

    Collects sequences from chai1_results.csv and/or native model outputs.
    """
    print("[Step3] Collecting sequences for analysis...")

    # Collect all unique candidates with sequences
    candidates: dict[str, dict[str, str]] = {}  # keyed by candidate_id

    # From chai1 results
    if chai1_csv and chai1_csv.exists():
        for row in _read_csv(chai1_csv):
            cid = row.get("candidate_id", "")
            seq = row.get("pred_sequence", "")
            if cid and seq:
                candidates[cid] = {
                    "candidate_id": cid,
                    "sample_id": row.get("sample_id", ""),
                    "model": row.get("model", ""),
                    "pred_sequence": seq,
                }

    # From native model outputs
    if native_roots:
        for root in native_roots:
            if not root.is_dir():
                continue
            for model_dir in sorted(root.iterdir()):
                if not model_dir.is_dir() or model_dir.name == "manifest.csv":
                    continue
                if models and model_dir.name not in models:
                    continue
                for sample_dir in sorted(model_dir.iterdir()):
                    if not sample_dir.is_dir():
                        continue
                    manifest = sample_dir / "candidate_manifest.csv"
                    if not manifest.exists():
                        continue
                    for row in _read_csv(manifest):
                        if row.get("status") != "ok":
                            continue
                        cid = row.get("candidate_id", "")
                        if not cid or cid in candidates:
                            continue

                        # Try to get sequence
                        seq = ""
                        seq_path = row.get("sequence_path", "")
                        if seq_path and Path(seq_path).exists():
                            try:
                                text = Path(seq_path).read_text(encoding="utf-8")
                                seq = "".join(
                                    l.strip() for l in text.splitlines()
                                    if not l.startswith(">")
                                )
                            except Exception:
                                pass
                        if not seq:
                            struct_path = row.get("structure_path", "")
                            if struct_path and Path(struct_path).exists():
                                parts = sample_dir.name.rsplit("_", 1)
                                hint = parts[-1] if len(parts) >= 2 else ""
                                seq = _extract_seq_from_structure(Path(struct_path), hint)

                        if seq:
                            candidates[cid] = {
                                "candidate_id": cid,
                                "sample_id": sample_dir.name,
                                "model": model_dir.name,
                                "pred_sequence": seq,
                            }

    print(f"[Step3] Analyzing {len(candidates)} sequences...")

    results: list[dict[str, str]] = []
    for i, (cid, cand) in enumerate(candidates.items(), 1):
        row = {k: "" for k in SEQUENCE_FIELDS}
        row["candidate_id"] = cid
        row["sample_id"] = cand["sample_id"]
        row["model"] = cand["model"]

        seq = cand["pred_sequence"].strip().upper()
        if seq:
            row["sequence_status"] = "ok"
            metrics = _analyze_sequence(seq)
            row.update(metrics)
        else:
            row["sequence_status"] = "no_sequence"

        results.append(row)

    _write_results(out_csv, results)
    print(f"[Step3] Done. {len(results)} results written to {out_csv}")
    return out_csv


def main() -> int:
    parser = argparse.ArgumentParser(description="Step 3: Sequence analysis")
    parser.add_argument(
        "--chai1-csv", type=Path,
        default=Path("outputs/evaluation/results/chai1_results.csv"),
    )
    parser.add_argument(
        "--native-root", action="append", type=Path, default=None,
    )
    parser.add_argument(
        "--index-csv", type=Path,
        default=Path("data/prepared/dataset_index_ready.csv"),
    )
    parser.add_argument(
        "--out-csv", type=Path,
        default=Path("outputs/evaluation/results/sequence_results.csv"),
    )
    parser.add_argument("--model", action="append", type=str, default=None)
    args = parser.parse_args()

    run_step3_sequence(
        chai1_csv=args.chai1_csv,
        native_roots=args.native_root,
        index_csv=args.index_csv,
        out_csv=args.out_csv,
        models=args.model,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
