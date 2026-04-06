#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path

from evaluation.schema import EXTERNAL_SEED_FIELDS


def _load_csv(path: Path) -> list[dict[str, str]]:
  if not path.exists():
    return []
  with path.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=EXTERNAL_SEED_FIELDS)
    writer.writeheader()
    writer.writerows(rows)


def _to_float(value: str) -> float | None:
  text = (value or "").strip()
  if not text:
    return None
  try:
    return float(text)
  except ValueError:
    return None


def _round_or_blank(value: float | None) -> str:
  if value is None:
    return ""
  return str(round(value, 6))


def _jaccard(a: set[str], b: set[str]) -> float:
  if not a and not b:
    return 1.0
  if not a or not b:
    return 0.0
  return len(a & b) / len(a | b)


def compute_robustness(in_csv: Path, out_csv: Path) -> Path:
  rows = _load_csv(in_csv)
  by_candidate: dict[str, list[dict[str, str]]] = {}
  for row in rows:
    cid = (row.get("candidate_id") or "").strip()
    if cid:
      by_candidate.setdefault(cid, []).append(row)

  out_rows: list[dict[str, str]] = []
  for cid, candidate_rows in by_candidate.items():
    pass_flags = [1.0 if (r.get("binding_plausibility_pass") or "") == "1" else 0.0 for r in candidate_rows]
    pass_rate = sum(pass_flags) / len(pass_flags) if pass_flags else 0.0

    pdockq2_values = [_to_float(r.get("pdockq2", "")) for r in candidate_rows]
    pdockq2_values = [v for v in pdockq2_values if v is not None]
    metric_std = statistics.pstdev(pdockq2_values) if len(pdockq2_values) > 1 else 0.0
    best_seed_vs_median_gap = None
    if pdockq2_values:
      med = statistics.median(pdockq2_values)
      best_seed_vs_median_gap = max(pdockq2_values) - med

    hotspot_sets = []
    for r in candidate_rows:
      labels = {x.strip() for x in (r.get("hotspot_string") or "").split(",") if x.strip()}
      hotspot_sets.append(labels)
    epitope_consistency = 1.0
    if len(hotspot_sets) > 1:
      pairwise = []
      for i in range(len(hotspot_sets)):
        for j in range(i + 1, len(hotspot_sets)):
          pairwise.append(_jaccard(hotspot_sets[i], hotspot_sets[j]))
      epitope_consistency = sum(pairwise) / len(pairwise) if pairwise else 1.0

    interface_residue_jaccard = epitope_consistency
    robust_pass = (pass_rate >= 0.67) and (metric_std <= 0.1)

    for row in candidate_rows:
      cur = {k: (row.get(k, "") or "") for k in EXTERNAL_SEED_FIELDS}
      cur["judge_seed_pass_rate"] = _round_or_blank(pass_rate)
      cur["judge_seed_metric_std"] = _round_or_blank(metric_std)
      cur["best_seed_vs_median_gap"] = _round_or_blank(best_seed_vs_median_gap)
      if not cur.get("secondary_judge_agreement"):
        cur["secondary_judge_agreement"] = ""
      cur["interface_residue_jaccard_across_seeds"] = _round_or_blank(interface_residue_jaccard)
      cur["epitope_consistency_across_seeds"] = _round_or_blank(epitope_consistency)
      cur["robust_pass"] = "1" if robust_pass else "0"
      out_rows.append(cur)

  _write_csv(out_csv, out_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="计算 Robustness 指标（基于 seed 聚合）")
  parser.add_argument(
    "--in-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="seed 级主表",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="输出 seed 级主表（可覆盖输入）",
  )
  args = parser.parse_args()
  out = compute_robustness(args.in_csv, args.out_csv)
  print(f"[OK] robustness: {out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
