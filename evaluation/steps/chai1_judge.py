#!/usr/bin/env python3
"""Step 1: Chai-1 structure prediction judge.

For each candidate, runs Chai-1 to predict the antibody-antigen complex
and extracts confidence metrics (iPTM, pTM, pLDDT).

Input:  Model output directories containing candidate_manifest.csv
Output: chai1_results.csv with continuous confidence metrics per candidate.

Supports resume: skips candidates whose result JSON already exists.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

from evaluation.structure_inference import load_reference_chain_sequence, protein_chain_sequences, parse_structure

# ---------------------------------------------------------------------------
# Output schema — intentionally minimal, all continuous metrics
# ---------------------------------------------------------------------------
CHAI1_FIELDS = [
    "candidate_id",
    "sample_id",
    "model",
    "pred_sequence",
    "chai1_status",
    "chai1_error",
    "chai1_iptm",
    "chai1_ptm",
    "chai1_plddt_binder",
    "chai1_structure_path",
    "chai1_confidence_path",
    "target_structure_path",
    "target_chain_id",
    "binder_chain_id",
]


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _write_results(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CHAI1_FIELDS)
        w.writeheader()
        w.writerows(rows)


def _extract_sequence_from_pdb(pdb_path: Path, chain_hint: str = "") -> str:
    """Extract sequence from PDB/CIF when no explicit sequence file exists."""
    try:
        structure = parse_structure(pdb_path)
        seqs = protein_chain_sequences(structure)
        if chain_hint and chain_hint in seqs:
            return seqs[chain_hint]
        if seqs:
            return min(seqs.values(), key=len)  # shortest = likely antibody
    except Exception:
        pass
    return ""


def _load_existing_results(result_json: Path) -> dict | None:
    if not result_json.exists():
        return None
    try:
        data = json.loads(result_json.read_text(encoding="utf-8"))
        if data.get("status") == "ok":
            return data
    except Exception:
        pass
    return None


def _collect_candidates_from_csv(
    candidates_csv: Path,
    index_csv: Path,
    models: list[str] | None = None,
) -> list[dict[str, str]]:
    """[antibench MVP] 从 benchmark_2 产出的 candidates_csv 读取候选列表。

    candidates_csv 字段：candidate_id, target_id, model_name, sequence,
                         structure_path, status, ...
    index_csv 字段：sample_id (= target_id), reference_complex_path,
                    antigen_chain, antibody_chain
    """
    # 加载 index 映射 target_id → 结构路径和链信息
    index_map: dict[str, dict[str, str]] = {}
    if index_csv.exists():
        for row in _read_csv(index_csv):
            sid = row.get("sample_id", "")
            if sid:
                index_map[sid] = row

    candidates = []
    for row in _read_csv(candidates_csv):
        status = row.get("status", "ok")
        if status == "failed":
            continue

        # source_type: native_reference 跳过（不需要再跑 Chai-1）
        if row.get("source_type", "generated") == "native_reference":
            continue

        target_id = row.get("target_id", "")
        model_name = row.get("model_name", "")
        if not target_id or not model_name:
            continue

        if models and model_name not in models:
            continue

        seq = row.get("sequence", "").strip().upper()
        struct_path = row.get("structure_path", "")
        candidate_id = row.get("candidate_id", "")

        # 从 index_csv 获取靶点结构路径和链 ID
        idx = index_map.get(target_id, {})
        target_structure_path = idx.get("reference_complex_path", "")
        antigen_chain = idx.get("antigen_chain", "")
        antibody_chain = idx.get("antibody_chain", "")

        # 若 target_structure_path 未在 index 中，尝试从旧字段回退
        if not target_structure_path:
            target_structure_path = idx.get("reference_structure_path", "")

        # 若序列为空，尝试从结构文件提取
        if not seq and struct_path and Path(struct_path).exists():
            ab_hint = antibody_chain or target_id.rsplit("_", 1)[-1]
            seq = _extract_sequence_from_pdb(Path(struct_path), ab_hint)

        candidates.append({
            "candidate_id": candidate_id,
            "sample_id": target_id,      # step1 内部用 sample_id
            "model": model_name,
            "pred_sequence": seq,
            "source_structure_path": struct_path,
            "target_structure_path": target_structure_path,
            "target_chain_id": antigen_chain,
            "binder_chain_id": antibody_chain,
        })

    return candidates


def _collect_candidates(
    native_roots: list[Path],
    index_csv: Path,
    models: list[str] | None = None,
) -> list[dict[str, str]]:
    """Scan native output directories for ok candidates and build a flat list."""
    # Load dataset index for target info
    index_map: dict[str, dict[str, str]] = {}
    if index_csv.exists():
        for row in _read_csv(index_csv):
            sid = row.get("sample_id", "")
            if sid:
                index_map[sid] = row

    candidates = []
    for root in native_roots:
        if not root.is_dir():
            continue
        # Structure: root/model_name/sample_id/candidate_manifest.csv
        for model_dir in sorted(root.iterdir()):
            if not model_dir.is_dir():
                continue
            model_name = model_dir.name
            if model_name == "manifest.csv":
                continue
            if models and model_name not in models:
                continue
            for sample_dir in sorted(model_dir.iterdir()):
                if not sample_dir.is_dir():
                    continue
                manifest = sample_dir / "candidate_manifest.csv"
                if not manifest.exists():
                    continue
                sample_id = sample_dir.name
                idx_row = index_map.get(sample_id, {})

                for row in _read_csv(manifest):
                    if row.get("status") != "ok":
                        continue
                    cid = row.get("candidate_id", "")
                    if not cid:
                        continue

                    # Get sequence: from manifest, from sequence file, or from structure
                    seq = ""
                    seq_path = row.get("sequence_path", "")
                    if seq_path and Path(seq_path).exists():
                        try:
                            text = Path(seq_path).read_text(encoding="utf-8")
                            seq = "".join(
                                line.strip()
                                for line in text.splitlines()
                                if not line.startswith(">")
                            )
                        except Exception:
                            pass

                    if not seq:
                        struct_path = row.get("structure_path", "")
                        if struct_path and Path(struct_path).exists():
                            # Antibody chain hint from sample_id (e.g. 7z1x_B_A → A is antibody)
                            parts = sample_id.rsplit("_", 1)
                            ab_hint = parts[-1] if len(parts) >= 2 else ""
                            seq = _extract_sequence_from_pdb(Path(struct_path), ab_hint)

                    # Target info
                    target_path = idx_row.get("reference_complex_path", "")
                    antigen_chain = idx_row.get("antigen_chain", "")
                    ab_chain = idx_row.get("antibody_chain", "")

                    candidates.append({
                        "candidate_id": cid,
                        "sample_id": sample_id,
                        "model": model_name,
                        "pred_sequence": seq,
                        "source_structure_path": row.get("structure_path", ""),
                        "target_structure_path": target_path,
                        "target_chain_id": antigen_chain,
                        "binder_chain_id": ab_chain,
                    })

    return candidates


def _run_chai1_single(
    candidate: dict[str, str],
    run_dir: Path,
    device: str = "cuda:0",
) -> dict[str, str]:
    """Run Chai-1 for a single candidate. Returns a result row."""
    cid = candidate["candidate_id"]
    result_row = {k: "" for k in CHAI1_FIELDS}
    result_row["candidate_id"] = cid
    result_row["sample_id"] = candidate["sample_id"]
    result_row["model"] = candidate["model"]
    result_row["pred_sequence"] = candidate["pred_sequence"]
    result_row["target_structure_path"] = candidate["target_structure_path"]
    result_row["target_chain_id"] = candidate["target_chain_id"]
    result_row["binder_chain_id"] = candidate["binder_chain_id"]

    seq = candidate["pred_sequence"]
    if not seq:
        result_row["chai1_status"] = "failed"
        result_row["chai1_error"] = "no_sequence"
        return result_row

    target_path = candidate["target_structure_path"]
    if not target_path or not Path(target_path).exists():
        result_row["chai1_status"] = "failed"
        result_row["chai1_error"] = "no_target_structure"
        return result_row

    # Prepare job directory
    job_dir = run_dir / cid.replace(":", "_").replace("/", "_")
    job_dir.mkdir(parents=True, exist_ok=True)

    # Check for existing result (resume support)
    result_json = job_dir / "result.json"
    existing = _load_existing_results(result_json)
    if existing:
        result_row["chai1_status"] = "ok"
        result_row["chai1_iptm"] = existing.get("iptm", "")
        result_row["chai1_ptm"] = existing.get("ptm", "")
        result_row["chai1_plddt_binder"] = existing.get("plddt_binder", "")
        result_row["chai1_structure_path"] = existing.get("structure_path", "")
        result_row["chai1_confidence_path"] = existing.get("confidence_path", "")
        return result_row

    # Prepare target structure (extract antigen chain)
    target_chain_hint = candidate["target_chain_id"]
    target_seq = load_reference_chain_sequence(
        Path(target_path), preferred_chain_id=target_chain_hint
    )
    if not target_seq:
        result_row["chai1_status"] = "failed"
        result_row["chai1_error"] = "cannot_extract_target_sequence"
        return result_row

    # Write FASTA
    fasta_path = job_dir / "input.fasta"
    fasta_body = (
        f">protein|name=target\n{target_seq}\n"
        f">protein|name=binder\n{seq}\n"
    )
    fasta_path.write_text(fasta_body, encoding="utf-8")

    # Run Chai-1 inference
    output_dir = job_dir / "chai_output"
    if output_dir.exists():
        import shutil
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    try:
        from chai_lab.chai1 import run_inference

        num_diffn_samples = int(os.environ.get("CHAI_NUM_DIFFN_SAMPLES", "5"))
        num_trunk_recycles = int(os.environ.get("CHAI_NUM_TRUNK_RECYCLES", "3"))
        num_diffn_timesteps = int(os.environ.get("CHAI_NUM_DIFFN_TIMESTEPS", "200"))
        use_msa_server = os.environ.get("CHAI_USE_MSA_SERVER", "").lower() in ("1", "true")
        low_memory = os.environ.get("CHAI_LOW_MEMORY", "1").lower() not in ("0", "false")

        candidates_out = run_inference(
            fasta_path,
            output_dir=output_dir,
            use_esm_embeddings=True,
            use_msa_server=use_msa_server,
            num_trunk_recycles=num_trunk_recycles,
            num_diffn_timesteps=num_diffn_timesteps,
            num_diffn_samples=num_diffn_samples,
            seed=42,
            device=device,
            low_memory=low_memory,
        )

        sorted_c = candidates_out.sorted()
        if not sorted_c.cif_paths:
            raise FileNotFoundError("chai_returned_no_structures")

        best_cif = sorted_c.cif_paths[0]
        best_rank = sorted_c.ranking_data[0]

        # Extract metrics
        iptm = ""
        try:
            iptm = str(round(float(best_rank.ptm_scores.interface_ptm.detach().cpu().item()), 6))
        except Exception:
            pass

        ptm = ""
        try:
            ptm = str(round(float(best_rank.ptm_scores.ptm.detach().cpu().item()), 6))
        except Exception:
            pass

        plddt_binder = ""
        try:
            per = best_rank.plddt_scores.per_chain_plddt.detach().cpu().float().numpy().ravel()
            if per.size > 0:
                val = float(per[-1])  # binder is last chain in FASTA
                if val > 1.5:
                    val = val / 100.0
                plddt_binder = str(round(val, 6))
        except Exception:
            pass

        # Save ranking data
        confidence_path = job_dir / "confidence.json"
        try:
            from chai_lab.ranking.rank import get_scores
            raw_scores = get_scores(best_rank)
            out_scores: dict = {}
            for k, v in raw_scores.items():
                out_scores[k] = v.tolist() if hasattr(v, "tolist") else v
            confidence_path.write_text(
                json.dumps(out_scores, indent=2, ensure_ascii=True) + "\n",
                encoding="utf-8",
            )
        except Exception:
            confidence_path.write_text("{}", encoding="utf-8")

        # Save result JSON for resume
        result_data = {
            "status": "ok",
            "iptm": iptm,
            "ptm": ptm,
            "plddt_binder": plddt_binder,
            "structure_path": str(best_cif.resolve()),
            "confidence_path": str(confidence_path.resolve()),
        }
        result_json.write_text(
            json.dumps(result_data, indent=2, ensure_ascii=True) + "\n",
            encoding="utf-8",
        )

        result_row["chai1_status"] = "ok"
        result_row["chai1_iptm"] = iptm
        result_row["chai1_ptm"] = ptm
        result_row["chai1_plddt_binder"] = plddt_binder
        result_row["chai1_structure_path"] = str(best_cif.resolve())
        result_row["chai1_confidence_path"] = str(confidence_path.resolve())

    except Exception as exc:
        result_row["chai1_status"] = "failed"
        result_row["chai1_error"] = f"{type(exc).__name__}:{exc}"
        # Save failure for resume
        result_json.write_text(
            json.dumps({"status": "failed", "error": str(exc)}) + "\n",
            encoding="utf-8",
        )

    return result_row


def run_step1_chai1(
    *,
    native_roots: list[Path],
    index_csv: Path,
    out_csv: Path,
    run_dir: Path,
    device: str = "cuda:0",
    models: list[str] | None = None,
    max_candidates: int = 0,
    # [antibench MVP] 若提供，优先从此 CSV 读取候选（benchmark_2 输出）
    candidates_csv: Path | None = None,
) -> Path:
    """Run Chai-1 evaluation on all candidates from native model outputs.

    Args:
        native_roots: Directories containing model outputs (model/sample/candidate_manifest.csv).
                      若 candidates_csv 已提供，此参数可为空列表。
        index_csv: Dataset index CSV with target structure paths and chain IDs.
        out_csv: Output CSV path for Chai-1 results.
        run_dir: Working directory for Chai-1 inference outputs.
        device: CUDA device string.
        models: Filter to specific model names, or None for all.
        max_candidates: Limit number of candidates (0 = no limit).
        candidates_csv: [antibench MVP] 若提供，从此 CSV 读取候选列表，
                        而非扫描 native_roots 目录。
    """
    run_dir.mkdir(parents=True, exist_ok=True)

    # ── 候选来源路由 ─────────────────────────────────────────────────────
    if candidates_csv is not None and candidates_csv.exists():
        print(f"[Step1] Reading candidates from {candidates_csv} ...")
        all_candidates = _collect_candidates_from_csv(candidates_csv, index_csv, models)
    else:
        print(f"[Step1] Collecting candidates from {len(native_roots)} root(s)...")
        all_candidates = _collect_candidates(native_roots, index_csv, models)

    if max_candidates > 0:
        all_candidates = all_candidates[:max_candidates]
    print(f"[Step1] Found {len(all_candidates)} candidates")

    results: list[dict[str, str]] = []

    # Load existing results for resume
    existing_results: dict[str, dict[str, str]] = {}
    if out_csv.exists():
        for row in _read_csv(out_csv):
            cid = row.get("candidate_id", "")
            if cid and row.get("chai1_status") == "ok":
                existing_results[cid] = row

    for i, cand in enumerate(all_candidates, 1):
        cid = cand["candidate_id"]

        # Resume: skip if already evaluated successfully
        if cid in existing_results:
            print(f"[Step1] ({i}/{len(all_candidates)}) {cid} — cached ✓")
            results.append(existing_results[cid])
            continue

        print(f"[Step1] ({i}/{len(all_candidates)}) {cid} — running Chai-1...")
        row = _run_chai1_single(cand, run_dir, device=device)
        results.append(row)

        # Write incrementally so progress is saved even on crash
        _write_results(out_csv, results)
        status = row["chai1_status"]
        iptm = row.get("chai1_iptm", "")
        suffix = f"iPTM={iptm}" if iptm else row.get("chai1_error", "")
        print(f"[Step1]   → {status}: {suffix}")

    _write_results(out_csv, results)
    print(f"[Step1] Done. {len(results)} results written to {out_csv}")
    return out_csv


def main() -> int:
    parser = argparse.ArgumentParser(description="Step 1: Run Chai-1 evaluation")
    parser.add_argument(
        "--native-root", action="append", type=Path, required=True,
        help="Model output directory (can specify multiple times)",
    )
    parser.add_argument(
        "--index-csv", type=Path,
        default=Path("data/prepared/dataset_index_ready.csv"),
    )
    parser.add_argument(
        "--out-csv", type=Path,
        default=Path("outputs/evaluation/results/chai1_results.csv"),
    )
    parser.add_argument(
        "--run-dir", type=Path,
        default=Path("outputs/evaluation/results/chai1_runs"),
    )
    parser.add_argument("--device", type=str, default="cuda:0")
    parser.add_argument("--model", action="append", type=str, default=None)
    parser.add_argument("--max-candidates", type=int, default=0)
    args = parser.parse_args()

    run_step1_chai1(
        native_roots=args.native_root,
        index_csv=args.index_csv,
        out_csv=args.out_csv,
        run_dir=args.run_dir,
        device=args.device,
        models=args.model,
        max_candidates=args.max_candidates,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
