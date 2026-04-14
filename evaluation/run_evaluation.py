#!/usr/bin/env python3
"""Unified evaluation pipeline — three-step architecture.

Usage:
  # Run all steps
  python -m evaluation.run_evaluation \
    --native-root outputs/final_boltzgen/boltzgen \
    --native-root outputs/native_predictions_cropped_formal/RFantibody \
    --out-dir outputs/evaluation/results

  # Run only Step 1 (Chai-1)
  python -m evaluation.run_evaluation --step 1 --native-root ...

  # Run only Steps 2+3 (no GPU needed if chai1 already done)
  python -m evaluation.run_evaluation --step 2 --step 3 --native-root ...

  # Run with PyRosetta interface analysis
  python -m evaluation.run_evaluation --use-pyrosetta --native-root ...

Steps:
  1. Chai-1 judge     → chai1_results.csv       (iPTM, pTM, pLDDT)
  2. Interface metrics → interface_results.csv   (BSA, pDockQ2, contacts, ddG)
  3. Sequence analysis → sequence_results.csv    (liabilities, CDR lengths, charge)
  → Aggregate         → eval_merged.csv         (all metrics in one table)
  4. Diversity+Epitope → diversity_epitope_results.csv (clustering, epitope coverage)
  → Final merge       → eval_final.csv          (all metrics + diversity + epitope)
"""
from __future__ import annotations

import argparse
from pathlib import Path


def run_evaluation(
    *,
    native_roots: list[Path],
    index_csv: Path,
    out_dir: Path,
    steps: list[int] | None = None,
    device: str = "cuda:0",
    models: list[str] | None = None,
    max_candidates: int = 0,
    use_pyrosetta: bool = False,
    epitope_dir: Path | None = None,
    # ── 新增：来自 antibench benchmark_2 的候选 CSV（可选） ─────────────────
    candidates_csv: Path | None = None,
) -> None:
    """Run the evaluation pipeline (all steps or selected ones).

    新增参数 candidates_csv：
        当通过 python -m antibench evaluate 调用时，
        由 benchmark_2_common_k 产出的候选 CSV 作为输入，
        而非扫描 native_roots 目录。
        若同时提供 native_roots，优先使用 candidates_csv。
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    chai1_csv = out_dir / "chai1_results.csv"
    interface_csv = out_dir / "interface_results.csv"
    sequence_csv = out_dir / "sequence_results.csv"
    merged_csv = out_dir / "eval_merged.csv"
    chai1_run_dir = out_dir / "chai1_runs"

    run_steps = set(steps) if steps else {1, 2, 3, 4}

    # Step 1: Chai-1
    if 1 in run_steps:
        from evaluation.steps.chai1_judge import run_step1_chai1

        run_step1_chai1(
            native_roots=native_roots,
            index_csv=index_csv,
            out_csv=chai1_csv,
            run_dir=chai1_run_dir,
            device=device,
            models=models,
            max_candidates=max_candidates,
            # [antibench MVP] 若提供 candidates_csv，Step 1 优先从中读取候选列表
            candidates_csv=candidates_csv,
        )

    # Step 2: Interface metrics
    if 2 in run_steps:
        from evaluation.steps.interface_metrics import run_step2_interface

        run_step2_interface(
            chai1_csv=chai1_csv if chai1_csv.exists() else None,
            native_roots=native_roots,
            index_csv=index_csv,
            out_csv=interface_csv,
            use_pyrosetta=use_pyrosetta,
        )

    # Step 3: Sequence analysis
    if 3 in run_steps:
        from evaluation.steps.sequence_metrics import run_step3_sequence

        run_step3_sequence(
            chai1_csv=chai1_csv if chai1_csv.exists() else None,
            native_roots=native_roots,
            index_csv=index_csv,
            out_csv=sequence_csv,
            models=models,
        )

    # Aggregate if any step ran
    if run_steps:
        from evaluation.aggregate import run_aggregate

        run_aggregate(
            chai1_csv=chai1_csv,
            interface_csv=interface_csv,
            sequence_csv=sequence_csv,
            out_csv=merged_csv,
        )

    # Step 4: Diversity + Epitope coverage
    diversity_csv = out_dir / "diversity_epitope_results.csv"
    final_csv = out_dir / "eval_final.csv"
    if 4 in run_steps:
        from evaluation.steps.diversity_epitope import run_step4_diversity_epitope

        # Auto-detect epitope dir if not specified
        epi = epitope_dir
        if epi is None:
            for candidate in [Path("data/epitopes"), Path("data/prepared/epitopes")]:
                if candidate.is_dir():
                    epi = candidate
                    break

        run_step4_diversity_epitope(
            merged_csv=merged_csv if merged_csv.exists() else out_dir / "eval_merged.csv",
            epitope_dir=epi,
            out_csv=diversity_csv,
        )

        # Merge diversity/epitope into final CSV
        from evaluation.aggregate import merge_extra_csv
        merge_extra_csv(
            base_csv=merged_csv,
            extra_csv=diversity_csv,
            out_csv=final_csv,
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Unified evaluation: Chai-1 + Interface + Sequence (three-step)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--native-root", action="append", type=Path, default=None,
        help="Model output directory containing sample subdirs (can repeat). "
             "At least one of --native-root or --candidates-csv is required.",
    )
    parser.add_argument(
        "--candidates-csv", type=Path, default=None,
        help="[antibench MVP] benchmark_2 产出的候选 CSV 路径，"
             "作为 evaluation pipeline 的统一输入替代 --native-root。",
    )
    parser.add_argument(
        "--index-csv", type=Path,
        default=Path("data/dataset_index_ready.csv"),
        help="Dataset index CSV",
    )
    parser.add_argument(
        "--out-dir", type=Path,
        default=Path("outputs/evaluation/results"),
        help="Output directory for all result CSVs",
    )
    parser.add_argument(
        "--step", action="append", type=int, default=None,
        help="Run only specific step(s): 1=Chai-1, 2=Interface, 3=Sequence",
    )
    parser.add_argument("--device", type=str, default="cuda:0")
    parser.add_argument(
        "--model", action="append", type=str, default=None,
        help="Filter to specific model names",
    )
    parser.add_argument("--max-candidates", type=int, default=0)
    parser.add_argument("--use-pyrosetta", action="store_true")
    parser.add_argument(
        "--epitope-dir", type=Path, default=None,
        help="Directory with {sample_id}.json epitope hotspot files",
    )
    args = parser.parse_args()

    run_evaluation(
        native_roots=args.native_root or [],
        index_csv=args.index_csv,
        out_dir=args.out_dir,
        steps=args.step,
        device=args.device,
        models=args.model,
        max_candidates=args.max_candidates,
        use_pyrosetta=args.use_pyrosetta,
        epitope_dir=args.epitope_dir,
        candidates_csv=args.candidates_csv,
    )
    print("[OK] Evaluation pipeline complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
