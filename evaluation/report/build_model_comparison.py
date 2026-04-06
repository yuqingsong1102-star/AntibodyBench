#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path

from Bio import Align

from evaluation.schema import EXTERNAL_MODEL_SUMMARY_FIELDS, EXTERNAL_SAMPLE_SUMMARY_FIELDS


def _load_csv(path: Path) -> list[dict[str, str]]:
  if not path.exists():
    return []
  with path.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
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


def _to_int(value: str) -> int | None:
  text = (value or "").strip()
  if not text:
    return None
  try:
    return int(float(text))
  except ValueError:
    return None


def _round_or_blank(value: float | None) -> str:
  if value is None:
    return ""
  return str(round(value, 6))


def _median(values: list[float]) -> float | None:
  return statistics.median(values) if values else None


def _sequence_identity(a: str, b: str, aligner: Align.PairwiseAligner) -> float:
  if not a or not b:
    return 0.0
  alignment = aligner.align(a, b)[0]
  return alignment.score / max(len(a), len(b))


def _accepted_diversity(rows: list[dict[str, str]]) -> float | None:
  accepted_sequences = [
    (row.get("pred_sequence") or "").strip()
    for row in rows
    if (_to_float(row.get("binding_pass_rate", "")) or 0.0) > 0.0 and (row.get("pred_sequence") or "").strip()
  ]
  if len(accepted_sequences) < 2:
    return None
  aligner = Align.PairwiseAligner()
  identities = []
  for i in range(len(accepted_sequences)):
    for j in range(i + 1, len(accepted_sequences)):
      identities.append(_sequence_identity(accepted_sequences[i], accepted_sequences[j], aligner))
  if not identities:
    return None
  return 1.0 - (sum(identities) / len(identities))


def build_model_comparison(
  *,
  sequence_csv: Path,
  seed_level_csv: Path,
  model_summary_csv: Path,
  sample_summary_csv: Path,
) -> tuple[Path, Path]:
  sequence_rows = _load_csv(sequence_csv)
  seed_rows = _load_csv(seed_level_csv)

  runtime_by_candidate: dict[str, float] = {}
  for row in seed_rows:
    cid = (row.get("candidate_id") or "").strip()
    if not cid or cid in runtime_by_candidate:
      continue
    duration = _to_float(row.get("source_duration_sec", ""))
    runtime_by_candidate[cid] = duration or 0.0

  by_model: dict[str, list[dict[str, str]]] = {}
  for row in sequence_rows:
    model = (row.get("model") or "").strip() or "unknown"
    by_model.setdefault(model, []).append(row)

  model_rows: list[dict[str, str]] = []
  for model in sorted(by_model):
    rows = by_model[model]
    by_sample: dict[str, list[dict[str, str]]] = {}
    for row in rows:
      sid = (row.get("sample_id") or "").strip() or "unknown"
      by_sample.setdefault(sid, []).append(row)

    def _pass_at_k(k: int) -> float | None:
      if not by_sample:
        return None
      hits = 0
      for sample_rows in by_sample.values():
        sample_rows = sorted(sample_rows, key=lambda x: _to_int(x.get("candidate_rank_in_model", "")) or 10**9)
        top_rows = sample_rows[:k]
        if any((_to_float(r.get("binding_pass_rate", "")) or 0.0) > 0.0 for r in top_rows):
          hits += 1
      return hits / len(by_sample)

    viability_rates = [_to_float(r.get("viability_pass_rate", "")) for r in rows]
    viability_rates = [v for v in viability_rates if v is not None]
    robust_count = sum(1 for r in rows if (r.get("robust_pass") or "") == "1")
    strong_binder_count = sum(1 for r in rows if (_to_float(r.get("affinity_proxy_pass_rate", "")) or 0.0) > 0.0)
    pdockq2_values = [_to_float(r.get("median_pdockq2", "")) for r in rows]
    pdockq2_values = [v for v in pdockq2_values if v is not None]
    iptm_values = [_to_float(r.get("median_external_iptm", "")) for r in rows]
    iptm_values = [v for v in iptm_values if v is not None]
    bsa_values = [_to_float(r.get("median_interface_bsa", "")) for r in rows]
    bsa_values = [v for v in bsa_values if v is not None]
    dg_values = [_to_float(r.get("median_interface_dG", "")) for r in rows]
    dg_values = [v for v in dg_values if v is not None]
    surf_values = [_to_float(r.get("median_surface_hydrophobicity", "")) for r in rows]
    surf_values = [v for v in surf_values if v is not None]

    accepted = [r for r in rows if (_to_float(r.get("binding_pass_rate", "")) or 0.0) > 0.0]
    total_runtime = sum(runtime_by_candidate.get((r.get("candidate_id") or ""), 0.0) for r in rows)
    runtime_per_accepted = (total_runtime / len(accepted)) if accepted else None

    summary = {field: "" for field in EXTERNAL_MODEL_SUMMARY_FIELDS}
    summary["model"] = model
    summary["n_candidates"] = str(len(rows))
    summary["viability_rate"] = _round_or_blank(_median(viability_rates))
    summary["bindability_pass_at_1"] = _round_or_blank(_pass_at_k(1))
    summary["bindability_pass_at_5"] = _round_or_blank(_pass_at_k(5))
    summary["bindability_pass_at_10"] = _round_or_blank(_pass_at_k(10))
    summary["strong_binder_proxy_rate"] = _round_or_blank(strong_binder_count / len(rows) if rows else None)
    summary["robust_pass_rate"] = _round_or_blank(robust_count / len(rows) if rows else None)
    summary["median_pdockq2"] = _round_or_blank(_median(pdockq2_values))
    summary["median_external_iptm"] = _round_or_blank(_median(iptm_values))
    summary["median_interface_bsa"] = _round_or_blank(_median(bsa_values))
    summary["median_interface_dG"] = _round_or_blank(_median(dg_values))
    summary["median_surface_hydrophobicity"] = _round_or_blank(_median(surf_values))
    summary["median_binding_plausibility_metrics"] = (
      f"pdockq2={summary['median_pdockq2']},iptm={summary['median_external_iptm']},bsa={summary['median_interface_bsa']}"
    )
    summary["median_affinity_proxy_metrics"] = (
      f"dG={summary['median_interface_dG']},surf_hydro={summary['median_surface_hydrophobicity']}"
    )
    summary["accepted_diversity"] = _round_or_blank(_accepted_diversity(rows))
    summary["runtime_per_accepted_binder"] = _round_or_blank(runtime_per_accepted)
    model_rows.append(summary)

  sample_rows: list[dict[str, str]] = []
  by_sample2: dict[str, list[dict[str, str]]] = {}
  for row in sequence_rows:
    sid = (row.get("sample_id") or "").strip() or "unknown"
    by_sample2.setdefault(sid, []).append(row)
  for sid in sorted(by_sample2):
    rows = by_sample2[sid]
    best_model = ""
    best_pdockq2 = None
    for row in rows:
      value = _to_float(row.get("best_pdockq2", ""))
      if value is None:
        continue
      if best_pdockq2 is None or value > best_pdockq2:
        best_pdockq2 = value
        best_model = (row.get("model") or "").strip()
    sample_rows.append(
      {
        "sample_id": sid,
        "n_models_attempted": str(len({(r.get('model') or '').strip() for r in rows if (r.get('model') or '').strip()})),
        "n_models_with_binding_hit": str(sum(1 for r in rows if (_to_float(r.get("binding_pass_rate", "")) or 0.0) > 0.0)),
        "best_model_by_pdockq": best_model,
        "best_pdockq": _round_or_blank(best_pdockq2),
      }
    )

  _write_csv(model_summary_csv, EXTERNAL_MODEL_SUMMARY_FIELDS, model_rows)
  _write_csv(sample_summary_csv, EXTERNAL_SAMPLE_SUMMARY_FIELDS, sample_rows)
  return model_summary_csv, sample_summary_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="构建 model/sample comparison 表")
  parser.add_argument(
    "--sequence-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_sequence_level.csv"),
    help="sequence 级长表",
  )
  parser.add_argument(
    "--seed-level-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="seed 级主表",
  )
  parser.add_argument(
    "--model-summary-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_model_summary.csv"),
    help="模型汇总表",
  )
  parser.add_argument(
    "--sample-summary-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_sample_summary.csv"),
    help="样本汇总表",
  )
  args = parser.parse_args()
  m, s = build_model_comparison(
    sequence_csv=args.sequence_csv,
    seed_level_csv=args.seed_level_csv,
    model_summary_csv=args.model_summary_csv,
    sample_summary_csv=args.sample_summary_csv,
  )
  print(f"[OK] 模型汇总: {m}")
  print(f"[OK] 样本汇总: {s}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
