#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path

from evaluation.schema import EXTERNAL_SEQUENCE_FIELDS


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


def _round_or_blank(value: float | None) -> str:
  if value is None:
    return ""
  return str(round(value, 6))


def _mean(values: list[float]) -> float | None:
  return (sum(values) / len(values)) if values else None


def _median(values: list[float]) -> float | None:
  return statistics.median(values) if values else None


def build_scorecards(seed_level_csv: Path, sequence_csv: Path) -> Path:
  seed_rows = _load_csv(seed_level_csv)
  by_candidate: dict[str, list[dict[str, str]]] = {}
  for row in seed_rows:
    cid = (row.get("candidate_id") or "").strip()
    if cid:
      by_candidate.setdefault(cid, []).append(row)

  out_rows: list[dict[str, str]] = []
  for cid in sorted(by_candidate):
    rows = by_candidate[cid]
    first = rows[0]

    viability_rates = [_to_float(r.get("viability_pass", "")) or 0.0 for r in rows]
    binding_rates = [_to_float(r.get("binding_plausibility_pass", "")) or 0.0 for r in rows]
    affinity_rates = [_to_float(r.get("affinity_proxy_pass", "")) or 0.0 for r in rows]
    developability_rates = [_to_float(r.get("developability_pass", "")) or 0.0 for r in rows]
    n_seed_success = sum(1 for r in rows if (r.get("external_prediction_ok") or "") == "1")

    pdockq2_values = [_to_float(r.get("pdockq2", "")) for r in rows]
    pdockq2_values = [v for v in pdockq2_values if v is not None]
    external_iptm_values = [_to_float(r.get("external_iptm", "")) for r in rows]
    external_iptm_values = [v for v in external_iptm_values if v is not None]
    external_ipae_values = [_to_float(r.get("external_ipae", "")) for r in rows]
    external_ipae_values = [v for v in external_ipae_values if v is not None]
    interface_residue_values = [_to_float(r.get("interface_residue_count", "")) for r in rows]
    interface_residue_values = [v for v in interface_residue_values if v is not None]
    antigen_contact_values = [_to_float(r.get("antigen_contact_residue_count", "")) for r in rows]
    antigen_contact_values = [v for v in antigen_contact_values if v is not None]
    contact_pair_values = [_to_float(r.get("interface_contact_pair_count", "")) for r in rows]
    contact_pair_values = [v for v in contact_pair_values if v is not None]
    bsa_values = [_to_float(r.get("interface_bsa", "")) for r in rows]
    bsa_values = [v for v in bsa_values if v is not None]
    dg_values = [_to_float(r.get("interface_dG", "")) for r in rows]
    dg_values = [v for v in dg_values if v is not None]
    hbond_values = [_to_float(r.get("interface_hbond_count", "")) for r in rows]
    hbond_values = [v for v in hbond_values if v is not None]
    clash_values = [_to_float(r.get("interface_clash_ratio", "")) for r in rows]
    clash_values = [v for v in clash_values if v is not None]
    hotspot_values = [_to_float(r.get("hotspot_contact_count", "")) for r in rows]
    hotspot_values = [v for v in hotspot_values if v is not None]
    binder_score_values = [_to_float(r.get("binder_score", "")) for r in rows]
    binder_score_values = [v for v in binder_score_values if v is not None]
    sap_values = [_to_float(r.get("sap_score", "")) for r in rows]
    sap_values = [v for v in sap_values if v is not None]
    surf_hydro_values = [_to_float(r.get("surface_hydrophobicity", "")) for r in rows]
    surf_hydro_values = [v for v in surf_hydro_values if v is not None]

    row = {field: "" for field in EXTERNAL_SEQUENCE_FIELDS}
    row["candidate_id"] = cid
    row["model"] = first.get("model", "")
    row["sample_id"] = first.get("sample_id", "")
    row["source_candidate_rank_in_model"] = first.get("source_candidate_rank_in_model", "")
    row["candidate_rank_in_model"] = first.get("candidate_rank_in_model", "")
    row["n_judge_seeds"] = str(len(rows))
    row["n_seed_success"] = str(n_seed_success)
    row["viability_pass_rate"] = _round_or_blank(_mean(viability_rates))
    row["binding_pass_rate"] = _round_or_blank(_mean(binding_rates))
    row["affinity_proxy_pass_rate"] = _round_or_blank(_mean(affinity_rates))
    row["developability_pass_rate"] = _round_or_blank(_mean(developability_rates))
    row["viability_gate_status"] = "pass" if (row["viability_pass_rate"] and float(row["viability_pass_rate"]) > 0.0) else "fail"
    row["robust_pass"] = "1" if (_mean(binding_rates) or 0.0) >= 0.67 else "0"
    row["best_pdockq2"] = _round_or_blank(max(pdockq2_values) if pdockq2_values else None)
    row["median_pdockq2"] = _round_or_blank(_median(pdockq2_values))
    row["median_external_iptm"] = _round_or_blank(_median(external_iptm_values))
    row["median_external_ipae"] = _round_or_blank(_median(external_ipae_values))
    row["median_interface_residue_count"] = _round_or_blank(_median(interface_residue_values))
    row["median_antigen_contact_residue_count"] = _round_or_blank(_median(antigen_contact_values))
    row["median_interface_contact_pair_count"] = _round_or_blank(_median(contact_pair_values))
    row["median_interface_bsa"] = _round_or_blank(_median(bsa_values))
    row["median_interface_dG"] = _round_or_blank(_median(dg_values))
    row["median_interface_hbond_count"] = _round_or_blank(_median(hbond_values))
    row["median_interface_clash_ratio"] = _round_or_blank(_median(clash_values))
    row["median_hotspot_contact_count"] = _round_or_blank(_median(hotspot_values))
    row["median_binder_score"] = _round_or_blank(_median(binder_score_values))
    row["median_sap_score"] = _round_or_blank(_median(sap_values))
    row["median_surface_hydrophobicity"] = _round_or_blank(_median(surf_hydro_values))
    row["accepted_set_diversity"] = ""
    row["sequence_entropy"] = first.get("sequence_entropy", "")
    row["any_basic_liability_flag"] = "1" if any((r.get("basic_liability_flag") or "") == "1" for r in rows) else "0"
    row["net_charge"] = first.get("net_charge", "")
    row["pred_seq_len"] = first.get("pred_seq_len", "")
    row["pred_sequence"] = first.get("pred_sequence", "")
    out_rows.append(row)

  _write_csv(sequence_csv, EXTERNAL_SEQUENCE_FIELDS, out_rows)
  return sequence_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="构建 sequence-level scorecards")
  parser.add_argument(
    "--seed-level-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="seed 级主表",
  )
  parser.add_argument(
    "--sequence-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_sequence_level.csv"),
    help="sequence 级长表",
  )
  args = parser.parse_args()
  out = build_scorecards(seed_level_csv=args.seed_level_csv, sequence_csv=args.sequence_csv)
  print(f"[OK] sequence-level: {out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
