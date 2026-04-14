#!/usr/bin/env python3
"""Step 2: Interface metrics from predicted complex structures.

Computes geometry-based interface quality metrics from PDB/CIF structures.
Uses BioPython by default. If PyRosetta is available and --use-pyrosetta is set,
computes additional Rosetta-level metrics (ddG, shape complementarity).

Input:  chai1_results.csv (from Step 1) — uses chai1_structure_path
        OR source model structures if Chai-1 wasn't run
Output: interface_results.csv with continuous interface metrics per candidate.
"""
from __future__ import annotations

import argparse
import copy
import csv
import math
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser
from Bio.PDB.SASA import ShrakeRupley

from evaluation.structure_inference import (
    load_reference_chain_sequence,
    protein_chain_sequences,
    resolve_structure_chain_ids,
)

# ---------------------------------------------------------------------------
# Output schema — continuous metrics only, no pass/fail
# ---------------------------------------------------------------------------
INTERFACE_FIELDS = [
    "candidate_id",
    "sample_id",
    "model",
    "interface_status",
    "interface_error",
    # Contact geometry
    "interface_residue_count",
    "antigen_contact_residue_count",
    "interface_contact_pair_count",
    "interface_clash_count",
    "interface_clash_ratio",
    # Buried surface area
    "interface_bsa",
    "binder_bsa_fraction",
    # pDockQ2 & pLDDT
    "pdockq2",
    "avg_interface_plddt",
    "binder_mean_plddt",
    # Surface properties
    "surface_hydrophobicity",
    "interface_fragmentation_index",
    # Hotspot contacts
    "hotspot_contact_count",
    # Rosetta metrics (empty if pyrosetta not used)
    "rosetta_ddG",
    "rosetta_sc",
    "rosetta_dSASA",
    "rosetta_hbonds_total",
    # Approximated energy terms (from BioPython)
    "approx_dG",
    "approx_hbond_count",
    "approx_saltbridge_count",
    # Source info
    "structure_source",  # "chai1" or "model_native"
    "structure_path",
]

CONTACT_CUTOFF = 5.0
CLASH_CUTOFF = 2.0
SURFACE_SASA_THRESHOLD = 15.0
HYDROPHOBIC_RES3 = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "PRO"}


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _write_results(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=INTERFACE_FIELDS)
        w.writeheader()
        w.writerows(rows)


def _r(value: float | None, digits: int = 6) -> str:
    return "" if value is None else str(round(value, digits))


def _parse_structure(path: Path):
    if path.suffix.lower() in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("s", str(path))


def _first_model(structure):
    return next(iter(structure))


def _get_chain(model, chain_id: str):
    for chain in model.get_chains():
        if chain.id == chain_id:
            return chain
    return None


def _heavy_atoms(chain) -> list:
    return [a for a in chain.get_atoms() if (getattr(a, "element", "") or "").strip().upper() != "H"]


def _extract_single_chain(structure, chain_id: str):
    s = copy.deepcopy(structure)
    m = _first_model(s)
    for c in list(m.get_chains()):
        if c.id != chain_id:
            m.detach_child(c.id)
    return s


def _sum_sasa(chain) -> float:
    return float(sum(float(getattr(a, "sasa", 0.0)) for a in _heavy_atoms(chain)))


def _compute_bsa(structure, target_cid: str, binder_cid: str) -> tuple[float | None, float | None]:
    try:
        sr = ShrakeRupley()
        full = copy.deepcopy(structure)
        sr.compute(full, level="A")
        fm = _first_model(full)
        t_bound = _sum_sasa(_get_chain(fm, target_cid))
        b_bound = _sum_sasa(_get_chain(fm, binder_cid))

        t_only = _extract_single_chain(structure, target_cid)
        b_only = _extract_single_chain(structure, binder_cid)
        sr.compute(t_only, level="A")
        sr.compute(b_only, level="A")
        t_free = _sum_sasa(_get_chain(_first_model(t_only), target_cid))
        b_free = _sum_sasa(_get_chain(_first_model(b_only), binder_cid))

        bsa = 0.5 * (max(0, t_free - t_bound) + max(0, b_free - b_bound))
        b_frac = ((b_free - b_bound) / b_free * 100) if b_free > 0 else None
        return bsa, b_frac
    except Exception:
        return None, None


def _compute_pdockq2(target_chain, binder_chain) -> tuple[float | None, float | None]:
    t_coords, t_plddt, b_coords, b_plddt = [], [], [], []
    for chain, coords, plddt in [(target_chain, t_coords, t_plddt), (binder_chain, b_coords, b_plddt)]:
        for res in chain.get_residues():
            if res.get_id()[0] != " ":
                continue
            atom = res["CB"] if "CB" in res else (res["CA"] if "CA" in res else None)
            if atom is None:
                continue
            coords.append(np.asarray(atom.get_coord(), dtype=float))
            plddt.append(float(atom.get_bfactor()))

    if not t_coords or not b_coords:
        return None, None

    dists = np.sqrt(np.sum(
        (np.asarray(t_coords)[:, None, :] - np.asarray(b_coords)[None, :, :]) ** 2,
        axis=-1,
    ))
    mask = dists <= 8.0
    n_contacts = int(np.count_nonzero(mask))
    if n_contacts == 0:
        return 0.0, 0.0

    t_iface = np.where(mask.any(axis=1))[0]
    b_iface = np.where(mask.any(axis=0))[0]
    plddt_vals = [t_plddt[i] for i in t_iface] + [b_plddt[i] for i in b_iface]
    if not plddt_vals or max(abs(v) for v in plddt_vals) < 1e-6:
        return None, None

    avg_plddt = float(np.mean(plddt_vals))
    x = avg_plddt * math.log10(n_contacts)
    pdockq2 = 0.724 / (1.0 + math.exp(-0.052 * (x - 152.611))) + 0.018
    return pdockq2, avg_plddt


def _compute_surface_hydrophobicity(structure, binder_cid: str) -> float | None:
    try:
        sr = ShrakeRupley()
        b_only = _extract_single_chain(structure, binder_cid)
        sr.compute(b_only, level="R")
        chain = _get_chain(_first_model(b_only), binder_cid)
        total, hydro = 0, 0
        for res in chain.get_residues():
            if res.get_id()[0] != " ":
                continue
            if float(getattr(res, "sasa", 0.0)) < SURFACE_SASA_THRESHOLD:
                continue
            total += 1
            if res.get_resname().strip().upper() in HYDROPHOBIC_RES3:
                hydro += 1
        return (hydro / total) if total > 0 else None
    except Exception:
        return None


def _mean_binder_plddt(chain) -> float | None:
    vals = []
    for res in chain.get_residues():
        if res.get_id()[0] != " ":
            continue
        atom = res["CA"] if "CA" in res else None
        if atom:
            vals.append(float(atom.get_bfactor()))
    if not vals or max(abs(v) for v in vals) < 1e-6:
        return None
    return float(np.mean(vals))


def _res_label(chain_id: str, residue) -> str:
    _, resseq, icode = residue.get_id()
    ic = "" if (icode or " ") == " " else icode
    return f"{chain_id}{resseq}{ic}"


def _compute_interface(
    structure_path: Path,
    target_seq: str,
    binder_seq: str,
    target_hint: str,
    binder_hint: str,
    hotspot_labels: set[str],
) -> dict[str, str]:
    """Compute all interface metrics from a single structure."""
    structure = _parse_structure(structure_path)
    chain_seqs = protein_chain_sequences(structure)

    t_cid, b_cid = resolve_structure_chain_ids(
        chain_seqs,
        target_sequence=target_seq,
        binder_sequence=binder_seq,
        target_chain_hint=target_hint,
        binder_chain_hint=binder_hint,
    )
    if not t_cid or not b_cid:
        raise ValueError(f"cannot_resolve_chains: target={t_cid} binder={b_cid}")

    model = _first_model(structure)
    t_chain = _get_chain(model, t_cid)
    b_chain = _get_chain(model, b_cid)
    if t_chain is None or b_chain is None:
        raise KeyError(f"missing_chain: target={t_cid} binder={b_cid}")

    # Contact analysis
    t_atoms = _heavy_atoms(t_chain)
    b_atoms = _heavy_atoms(b_chain)
    ns = NeighborSearch(t_atoms + b_atoms)
    t_ids = {id(a) for a in t_atoms}

    contact_pairs: set[tuple[str, str]] = set()
    t_contact_res: set[str] = set()
    b_contact_res: set[str] = set()
    clash_count = 0

    for ba in b_atoms:
        for ta in ns.search(ba.get_coord(), CONTACT_CUTOFF, level="A"):
            if id(ta) not in t_ids:
                continue
            tl = _res_label(t_cid, ta.get_parent())
            bl = _res_label(b_cid, ba.get_parent())
            contact_pairs.add((tl, bl))
            t_contact_res.add(tl)
            b_contact_res.add(bl)
            if float(ta - ba) < CLASH_CUTOFF:
                clash_count += 1

    n_pairs = len(contact_pairs)
    b_interface_count = len(b_contact_res)
    t_contact_count = len(t_contact_res)
    clash_ratio = (clash_count / max(1, n_pairs)) if n_pairs else None

    # Fragmentation
    b_resseqs = sorted({
        int(res.get_id()[1])
        for res in b_chain.get_residues()
        if res.get_id()[0] == " "
    })
    segs, prev = 0, None
    for rs in b_resseqs:
        if prev is None or rs > prev + 1:
            segs += 1
        prev = rs
    frag_idx = (segs / b_interface_count) if b_interface_count else None

    # BSA
    bsa, b_frac = _compute_bsa(structure, t_cid, b_cid)

    # pDockQ2
    pdockq2, avg_iface_plddt = _compute_pdockq2(t_chain, b_chain)

    # Surface hydrophobicity
    surf_hydro = _compute_surface_hydrophobicity(structure, b_cid)

    # Binder pLDDT
    binder_plddt = _mean_binder_plddt(b_chain)

    # Hotspot
    hotspot_count = len(t_contact_res & hotspot_labels) if hotspot_labels else 0

    # Approximate energy (from contact geometry)
    approx_hbond = max(0, int(round(n_pairs * 0.08)))
    approx_salt = max(0, int(round(n_pairs * 0.02)))
    approx_dG = None
    if bsa is not None:
        clash_pen = (clash_ratio or 0) * 50.0
        approx_dG = -0.018 * bsa - 0.25 * approx_hbond - 0.4 * approx_salt + clash_pen

    return {
        "interface_residue_count": str(b_interface_count),
        "antigen_contact_residue_count": str(t_contact_count),
        "interface_contact_pair_count": str(n_pairs),
        "interface_clash_count": str(clash_count),
        "interface_clash_ratio": _r(clash_ratio),
        "interface_bsa": _r(bsa),
        "binder_bsa_fraction": _r(b_frac),
        "pdockq2": _r(pdockq2),
        "avg_interface_plddt": _r(avg_iface_plddt),
        "binder_mean_plddt": _r(binder_plddt),
        "surface_hydrophobicity": _r(surf_hydro),
        "interface_fragmentation_index": _r(frag_idx),
        "hotspot_contact_count": str(hotspot_count),
        "approx_dG": _r(approx_dG),
        "approx_hbond_count": str(approx_hbond),
        "approx_saltbridge_count": str(approx_salt),
    }


def _run_pyrosetta_interface(pdb_path: Path, target_cid: str, binder_cid: str) -> dict[str, str]:
    """Run PyRosetta InterfaceAnalyzer if available. Returns empty strings if not."""
    try:
        import pyrosetta
        from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

        pyrosetta.init("-mute all", silent=True)
        pose = pyrosetta.pose_from_file(str(pdb_path))

        # Determine interface string (chains separated by _)
        interface = f"{target_cid}_{binder_cid}"
        iam = InterfaceAnalyzerMover(interface)
        iam.set_pack_separated(True)
        iam.set_compute_interface_sc(True)
        iam.apply(pose)

        return {
            "rosetta_ddG": _r(iam.get_interface_dG()),
            "rosetta_sc": _r(iam.get_interface_sc()),
            "rosetta_dSASA": _r(iam.get_interface_delta_sasa()),
            "rosetta_hbonds_total": str(iam.get_all_per_res_data().size()),
        }
    except ImportError:
        return {"rosetta_ddG": "", "rosetta_sc": "", "rosetta_dSASA": "", "rosetta_hbonds_total": ""}
    except Exception as exc:
        return {"rosetta_ddG": "", "rosetta_sc": "", "rosetta_dSASA": "", "rosetta_hbonds_total": ""}


def _analyze_one(
    candidate_id: str,
    sample_id: str,
    model_name: str,
    structure_path: Path,
    structure_source: str,
    target_path: str,
    target_hint: str,
    binder_hint: str,
    binder_seq: str,
    hotspot_string: str,
    use_pyrosetta: bool,
) -> dict[str, str]:
    """Analyze a single candidate structure."""
    row = {k: "" for k in INTERFACE_FIELDS}
    row["candidate_id"] = candidate_id
    row["sample_id"] = sample_id
    row["model"] = model_name
    row["structure_source"] = structure_source
    row["structure_path"] = str(structure_path)

    if not structure_path.exists():
        row["interface_status"] = "failed"
        row["interface_error"] = "structure_not_found"
        return row

    # Parse hotspot labels
    hotspot_labels = set()
    if hotspot_string:
        hotspot_labels = {s.strip() for s in hotspot_string.split(",") if s.strip()}

    # Get target sequence
    target_seq = ""
    if target_path and Path(target_path).exists():
        target_seq = load_reference_chain_sequence(Path(target_path), preferred_chain_id=target_hint)

    try:
        metrics = _compute_interface(
            structure_path, target_seq, binder_seq,
            target_hint, binder_hint, hotspot_labels,
        )
        row.update(metrics)
        row["interface_status"] = "ok"
    except Exception as exc:
        row["interface_status"] = "failed"
        row["interface_error"] = f"{type(exc).__name__}:{exc}"
        return row

    # Optional PyRosetta
    if use_pyrosetta and structure_path.suffix.lower() in {".pdb"}:
        rosetta_metrics = _run_pyrosetta_interface(structure_path, target_hint, binder_hint)
        row.update(rosetta_metrics)

    return row


def run_step2_interface(
    *,
    chai1_csv: Path | None = None,
    native_roots: list[Path] | None = None,
    index_csv: Path,
    out_csv: Path,
    use_pyrosetta: bool = False,
    prefer_chai1_structure: bool = True,
) -> Path:
    """Run interface metric computation.

    Can use Chai-1 predicted structures (from Step 1) or model native structures.
    """
    print(f"[Step2] Computing interface metrics...")

    # Load dataset index
    index_map: dict[str, dict[str, str]] = {}
    if index_csv.exists():
        for row in _read_csv(index_csv):
            sid = row.get("sample_id", "")
            if sid:
                index_map[sid] = row

    # Build work items: (candidate_id, structure_path, source, metadata)
    work: list[dict] = []

    if chai1_csv and chai1_csv.exists():
        for row in _read_csv(chai1_csv):
            if row.get("chai1_status") != "ok":
                continue
            struct = row.get("chai1_structure_path", "")
            if not struct or not Path(struct).exists():
                continue
            sid = row.get("sample_id", "")
            idx = index_map.get(sid, {})
            work.append({
                "candidate_id": row["candidate_id"],
                "sample_id": sid,
                "model": row.get("model", ""),
                "structure_path": Path(struct),
                "structure_source": "chai1",
                "target_path": row.get("target_structure_path", idx.get("reference_complex_path", "")),
                "target_hint": row.get("target_chain_id", idx.get("antigen_chain", "")),
                "binder_hint": row.get("binder_chain_id", idx.get("antibody_chain", "")),
                "binder_seq": row.get("pred_sequence", ""),
                "hotspot": idx.get("hotspot_string", ""),
            })

    # Also collect from native model structures if requested
    if native_roots and not prefer_chai1_structure:
        seen = {w["candidate_id"] for w in work}
        for root in native_roots:
            if not root.is_dir():
                continue
            for model_dir in sorted(root.iterdir()):
                if not model_dir.is_dir() or model_dir.name == "manifest.csv":
                    continue
                for sample_dir in sorted(model_dir.iterdir()):
                    if not sample_dir.is_dir():
                        continue
                    manifest = sample_dir / "candidate_manifest.csv"
                    if not manifest.exists():
                        continue
                    sid = sample_dir.name
                    idx = index_map.get(sid, {})
                    for row in _read_csv(manifest):
                        if row.get("status") != "ok":
                            continue
                        cid = row.get("candidate_id", "")
                        if not cid or cid in seen:
                            continue
                        struct = row.get("structure_path", "")
                        if not struct or not Path(struct).exists():
                            continue
                        work.append({
                            "candidate_id": cid,
                            "sample_id": sid,
                            "model": model_dir.name,
                            "structure_path": Path(struct),
                            "structure_source": "model_native",
                            "target_path": idx.get("reference_complex_path", ""),
                            "target_hint": idx.get("antigen_chain", ""),
                            "binder_hint": idx.get("antibody_chain", ""),
                            "binder_seq": "",
                            "hotspot": idx.get("hotspot_string", ""),
                        })
                        seen.add(cid)

    print(f"[Step2] {len(work)} structures to analyze")

    results: list[dict[str, str]] = []
    for i, item in enumerate(work, 1):
        cid = item["candidate_id"]
        print(f"[Step2] ({i}/{len(work)}) {cid}...")
        row = _analyze_one(
            candidate_id=cid,
            sample_id=item["sample_id"],
            model_name=item["model"],
            structure_path=item["structure_path"],
            structure_source=item["structure_source"],
            target_path=item["target_path"],
            target_hint=item["target_hint"],
            binder_hint=item["binder_hint"],
            binder_seq=item["binder_seq"],
            hotspot_string=item["hotspot"],
            use_pyrosetta=use_pyrosetta,
        )
        results.append(row)
        bsa = row.get("interface_bsa", "")
        pq = row.get("pdockq2", "")
        print(f"[Step2]   → {row['interface_status']}: BSA={bsa} pDockQ2={pq}")

    _write_results(out_csv, results)
    print(f"[Step2] Done. {len(results)} results written to {out_csv}")
    return out_csv


def main() -> int:
    parser = argparse.ArgumentParser(description="Step 2: Compute interface metrics")
    parser.add_argument(
        "--chai1-csv", type=Path,
        default=Path("outputs/evaluation/results/chai1_results.csv"),
        help="Chai-1 results from Step 1",
    )
    parser.add_argument(
        "--native-root", action="append", type=Path, default=None,
        help="Model output directory (optional, for native structures)",
    )
    parser.add_argument(
        "--index-csv", type=Path,
        default=Path("data/prepared/dataset_index_ready.csv"),
    )
    parser.add_argument(
        "--out-csv", type=Path,
        default=Path("outputs/evaluation/results/interface_results.csv"),
    )
    parser.add_argument("--use-pyrosetta", action="store_true")
    parser.add_argument(
        "--prefer-native", action="store_true",
        help="Use model native structures instead of Chai-1 structures",
    )
    args = parser.parse_args()

    run_step2_interface(
        chai1_csv=args.chai1_csv,
        native_roots=args.native_root,
        index_csv=args.index_csv,
        out_csv=args.out_csv,
        use_pyrosetta=args.use_pyrosetta,
        prefer_chai1_structure=not args.prefer_native,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
