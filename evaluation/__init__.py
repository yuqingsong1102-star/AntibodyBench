"""Evaluation package for AntibodyBench.

Layered layout:
- ingest/: build unified long tables from native outputs
- metrics/: compute per-run metrics on the long table
- report/: build summaries, figures, and markdown reports
- pipeline/: orchestrate the full evaluation workflow
"""

