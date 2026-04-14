"""Evaluation steps: three independent metric computation modules.

Step 1 (chai1_judge)      — Chai-1 structure prediction → iPTM, pTM, pLDDT
Step 2 (interface_metrics) — Interface geometry from predicted structures → BSA, contacts, pDockQ2, ddG (optional Rosetta)
Step 3 (sequence_metrics)  — Pure sequence analysis → liabilities, charge, entropy, CDR lengths
"""
