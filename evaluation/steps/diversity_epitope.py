#!/usr/bin/env python3
"""Step 4: Sequence diversity (clustering) and epitope coverage analysis.

Diversity:
  Uses MMseqs2 to cluster candidate sequences per (sample_id, model) at 90% identity.
  Reports cluster count, largest cluster fraction, and singleton fraction.

Epitope coverage:
  From Chai-1 predicted complex structures, extracts CDR-antigen contacts (<5 Å)
  and compares against known epitope hotspot residues.
  Reports coverage = |contacted hotspots| / |total hotspots|.

Input:  eval_merged.csv (from Steps 1-3)
Output: diversity_epitope_results.csv — one row per candidate
"""
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Output schema
# ---------------------------------------------------------------------------
DIVERSITY_EPITOPE_FIELDS = [
    "candidate_id",
    "sample_id",
    "model",
    # Diversity (per sample_id × model group)
    "cluster_id",          # which cluster this candidate belongs to
    "cluster_size",        # size of the cluster this candidate is in
    "group_cluster_count", # total clusters for this (sample, model)
    "group_total_seqs",    # total sequences in this group
    "group_largest_cluster_frac",  # largest cluster / total
    "group_singleton_frac",        # singletons / total
    # Epitope coverage (per candidate, requires Chai-1 structure)
    "epitope_contact_count",   # CDR residues contacting antigen hotspot
    "epitope_total_hotspots",  # total hotspot residues
    "epitope_coverage",        # contact_count / total_hotspots
]

CONTACT_CUTOFF = 5.0  # Angstroms


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _write_results(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=DIVERSITY_EPITOPE_FIELDS)
        w.writeheader()
        w.writerows(rows)


# ── Diversity via MMseqs2 ────────────────────────────────────

def _find_mmseqs() -> str | None:
    """Find mmseqs2 binary."""
    for p in [
        shutil.which("mmseqs"),
        os.path.expanduser("~/.local/bin/mmseqs"),
    ]:
        if p and os.path.isfile(p):
            return p
    return None


def _cluster_sequences(
    sequences: dict[str, str],  # candidate_id → sequence
    identity: float = 0.9,
    coverage: float = 0.8,
) -> dict[str, int] | None:
    """Cluster sequences with MMseqs2. Returns {candidate_id: cluster_id} or None on failure."""
    mmseqs = _find_mmseqs()
    if not mmseqs or len(sequences) < 2:
        return None

    with tempfile.TemporaryDirectory(prefix="mmseqs_") as tmpdir:
        fasta = Path(tmpdir) / "input.fasta"
        with fasta.open("w") as f:
            for cid, seq in sequences.items():
                f.write(f">{cid}\n{seq}\n")

        db = Path(tmpdir) / "db"
        clu = Path(tmpdir) / "clu"
        tsv = Path(tmpdir) / "clu.tsv"
        tmp = Path(tmpdir) / "tmp"
        tmp.mkdir()

        try:
            subprocess.run(
                [mmseqs, "createdb", str(fasta), str(db)],
                capture_output=True, check=True, timeout=120,
            )
            subprocess.run(
                [mmseqs, "cluster", str(db), str(clu), str(tmp),
                 "--min-seq-id", str(identity),
                 "-c", str(coverage),
                 "--cov-mode", "0"],
                capture_output=True, check=True, timeout=120,
            )
            subprocess.run(
                [mmseqs, "createtsv", str(db), str(db), str(clu), str(tsv)],
                capture_output=True, check=True, timeout=120,
            )

            # Parse TSV: representative\tmember
            cluster_map: dict[str, list[str]] = defaultdict(list)
            for line in tsv.read_text().strip().splitlines():
                parts = line.split("\t")
                if len(parts) >= 2:
                    rep, member = parts[0], parts[1]
                    cluster_map[rep].append(member)

            # Assign cluster IDs
            result: dict[str, int] = {}
            for idx, (rep, members) in enumerate(sorted(cluster_map.items()), 1):
                for m in members:
                    result[m] = idx
            return result

        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            return None


def _compute_group_diversity(
    cluster_assignments: dict[str, int],
    total_seqs: int,
) -> tuple[int, float, float]:
    """Returns (cluster_count, largest_cluster_frac, singleton_frac)."""
    if not cluster_assignments:
        return 0, 0.0, 0.0

    from collections import Counter
    sizes = Counter(cluster_assignments.values())
    n_clusters = len(sizes)
    largest = max(sizes.values())
    singletons = sum(1 for s in sizes.values() if s == 1)

    return (
        n_clusters,
        largest / total_seqs if total_seqs > 0 else 0.0,
        singletons / n_clusters if n_clusters > 0 else 0.0,
    )


# ── Epitope coverage ────────────────────────────────────────

def _load_epitope_hotspots(
    epitope_dir: Path,
    sample_id: str,
) -> set[int] | None:
    """Load hotspot residue indices from epitope JSON file."""
    json_path = epitope_dir / f"{sample_id}.json"
    if not json_path.exists():
        return None
    try:
        data = json.loads(json_path.read_text(encoding="utf-8"))
        hotspots = set()
        if isinstance(data, dict):
            # Format: {"hotspot_residues": [1,2,3,...]} or {"residues": [...]}
            for key in ("hotspot_residues", "residues", "epitope_residues"):
                if key in data and isinstance(data[key], list):
                    hotspots.update(int(r) for r in data[key])
                    break
            # Alternative: flat list of residue indices
            if not hotspots and "hotspot_indices" in data:
                hotspots.update(int(r) for r in data["hotspot_indices"])
        elif isinstance(data, list):
            hotspots.update(int(r) for r in data)
        return hotspots if hotspots else None
    except Exception:
        return None


def _compute_epitope_coverage(
    structure_path: Path,
    antigen_chain: str,
    binder_chain: str,
    hotspot_residues: set[int],
) -> tuple[int, float] | None:
    """Compute how many hotspot residues are contacted by the binder.

    Returns (contact_count, coverage) or None on failure.
    """
    try:
        from Bio.PDB import MMCIFParser, PDBParser

        suffix = structure_path.suffix.lower()
        if suffix in (".cif", ".mmcif"):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        structure = parser.get_structure("s", str(structure_path))
        model = next(iter(structure))

        ag_chain = None
        bd_chain = None
        for c in model.get_chains():
            if c.id == antigen_chain:
                ag_chain = c
            elif c.id == binder_chain:
                bd_chain = c

        # If exact chain IDs don't match, try by order (Chai-1 relabels chains)
        chains = list(model.get_chains())
        if ag_chain is None and len(chains) >= 1:
            ag_chain = chains[0]  # target is first in FASTA
        if bd_chain is None and len(chains) >= 2:
            bd_chain = chains[1]  # binder is second

        if ag_chain is None or bd_chain is None:
            return None

        # Get binder heavy atom coords
        binder_atoms = []
        for atom in bd_chain.get_atoms():
            if atom.element.strip().upper() != "H":
                binder_atoms.append(np.asarray(atom.get_coord(), dtype=float))

        if not binder_atoms:
            return None
        binder_coords = np.array(binder_atoms)

        # Check which antigen hotspot residues are contacted
        contacted = 0
        for res in ag_chain.get_residues():
            if res.get_id()[0] != " ":
                continue
            res_idx = res.get_id()[1]
            if res_idx not in hotspot_residues:
                continue
            # Check if any binder atom is within cutoff
            for atom in res.get_atoms():
                if atom.element.strip().upper() == "H":
                    continue
                coord = np.asarray(atom.get_coord(), dtype=float)
                dists = np.sqrt(np.sum((binder_coords - coord) ** 2, axis=1))
                if np.min(dists) <= CONTACT_CUTOFF:
                    contacted += 1
                    break

        coverage = contacted / len(hotspot_residues) if hotspot_residues else 0.0
        return contacted, coverage

    except Exception:
        return None


# ── Main entry ──────────────────────────────────────────────

def run_step4_diversity_epitope(
    *,
    merged_csv: Path,
    epitope_dir: Path | None = None,
    out_csv: Path,
) -> Path:
    """Run diversity clustering and epitope coverage analysis.

    Args:
        merged_csv: eval_merged.csv from aggregate step.
        epitope_dir: Directory containing {sample_id}.json epitope files.
        out_csv: Output CSV path.
    """
    if not merged_csv.exists():
        print("[Step4] No merged CSV found, skipping diversity/epitope analysis.")
        return out_csv

    rows = _read_csv(merged_csv)
    if not rows:
        print("[Step4] No candidates in merged CSV.")
        return out_csv

    print(f"[Step4] Analyzing {len(rows)} candidates...")

    # ── Phase 1: Diversity clustering ──
    # Group sequences by (sample_id, model)
    groups: dict[tuple[str, str], dict[str, str]] = defaultdict(dict)
    for r in rows:
        sid = r.get("sample_id", "")
        model = r.get("model", "")
        cid = r.get("candidate_id", "")
        seq = r.get("pred_sequence", "")
        if sid and model and cid and seq:
            groups[(sid, model)][cid] = seq

    # Cluster each group
    print(f"[Step4] Clustering {len(groups)} (sample, model) groups...")
    all_clusters: dict[str, dict] = {}  # cid → cluster info
    for (sid, model), seqs in sorted(groups.items()):
        cluster_map = _cluster_sequences(seqs)
        n_clusters, largest_frac, singleton_frac = (0, 0.0, 0.0)
        if cluster_map:
            n_clusters, largest_frac, singleton_frac = _compute_group_diversity(
                cluster_map, len(seqs)
            )
        for cid, seq in seqs.items():
            clu_id = cluster_map.get(cid, 0) if cluster_map else 0
            # Count cluster size
            clu_size = 0
            if cluster_map:
                clu_size = sum(1 for v in cluster_map.values() if v == clu_id)
            all_clusters[cid] = {
                "cluster_id": str(clu_id),
                "cluster_size": str(clu_size),
                "group_cluster_count": str(n_clusters),
                "group_total_seqs": str(len(seqs)),
                "group_largest_cluster_frac": f"{largest_frac:.4f}",
                "group_singleton_frac": f"{singleton_frac:.4f}",
            }
        status = f"{n_clusters} clusters" if cluster_map else "mmseqs2 unavailable"
        print(f"[Step4]   {sid}/{model}: {len(seqs)} seqs → {status}")

    # ── Phase 2: Epitope coverage ──
    epi_results: dict[str, dict] = {}
    if epitope_dir and epitope_dir.is_dir():
        print(f"[Step4] Computing epitope coverage from {epitope_dir}...")
        for r in rows:
            cid = r.get("candidate_id", "")
            sid = r.get("sample_id", "")
            struct_path = r.get("chai1_structure_path", "") or r.get("structure_path", "")
            ag_chain = r.get("target_chain_id", "")
            bd_chain = r.get("binder_chain_id", "")

            if not (cid and struct_path and Path(struct_path).exists()):
                continue

            hotspots = _load_epitope_hotspots(epitope_dir, sid)
            if hotspots is None:
                continue

            result = _compute_epitope_coverage(
                Path(struct_path), ag_chain, bd_chain, hotspots,
            )
            if result is not None:
                contacted, coverage = result
                epi_results[cid] = {
                    "epitope_contact_count": str(contacted),
                    "epitope_total_hotspots": str(len(hotspots)),
                    "epitope_coverage": f"{coverage:.4f}",
                }
        print(f"[Step4]   Computed epitope coverage for {len(epi_results)} candidates")
    else:
        print("[Step4] No epitope directory provided, skipping epitope coverage.")

    # ── Build output ──
    results: list[dict[str, str]] = []
    for r in rows:
        cid = r.get("candidate_id", "")
        out_row = {k: "" for k in DIVERSITY_EPITOPE_FIELDS}
        out_row["candidate_id"] = cid
        out_row["sample_id"] = r.get("sample_id", "")
        out_row["model"] = r.get("model", "")

        # Diversity
        if cid in all_clusters:
            out_row.update(all_clusters[cid])

        # Epitope
        if cid in epi_results:
            out_row.update(epi_results[cid])

        results.append(out_row)

    _write_results(out_csv, results)
    print(f"[Step4] Done. {len(results)} results written to {out_csv}")
    return out_csv
