"""Evaluation package for AntibodyBench — three-step architecture.

Steps:
  1. chai1_judge      — Chai-1 structure prediction (iPTM, pTM, pLDDT)
  2. interface_metrics — Interface geometry from predicted structures (BSA, pDockQ2, contacts)
  3. sequence_metrics  — Pure sequence analysis (liabilities, CDR lengths, charge)

Modules:
  run_evaluation.py     — Unified CLI entry point (--step 1/2/3)
  aggregate.py          — Merge step CSVs into eval_merged.csv
  schema.py             — CSV field definitions
  structure_inference.py — BioPython structure parsing utilities

See also: analysis/analyze_results.py for visualization and figures.
"""

