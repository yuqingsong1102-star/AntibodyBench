#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path

from schema import EVALUATION_LONG_FIELDS, SUMMARY_BY_MODEL_FIELDS, SUMMARY_BY_SAMPLE_FIELDS


AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def _to_float(v: str) -> float | None:
  s = (v or "").strip()
  if not s:
    return None
  try:
    return float(s)
  except ValueError:
    return None


def _to_int(v: str) -> int | None:
  s = (v or "").strip()
  if not s:
    return None
  try:
    return int(float(s))
  except ValueError:
    return None


def _percentile(values: list[float], q: float) -> float | None:
  if not values:
    return None
  vals = sorted(values)
  if len(vals) == 1:
    return vals[0]
  pos = (len(vals) - 1) * q
  low = int(pos)
  high = min(low + 1, len(vals) - 1)
  if low == high:
    return vals[low]
  frac = pos - low
  return vals[low] * (1.0 - frac) + vals[high] * frac


def _round_or_blank(x: float | None) -> str:
  if x is None:
    return ""
  return str(round(x, 6))


def _compute_seq_fields(row: dict[str, str]) -> dict[str, str]:
  seq = (row.get("pred_sequence") or "").strip().upper()
  if not seq:
    row["pred_seq_len"] = ""
    row["aa_valid_ratio"] = ""
    row["cys_count"] = ""
    row["motif_h3_len_in_8_20"] = ""
    row["has_sequence"] = "0"
    return row

  n = len(seq)
  valid_n = sum(1 for c in seq if c in AA20)
  cys_count = seq.count("C")
  row["pred_seq_len"] = str(n)
  row["aa_valid_ratio"] = _round_or_blank(valid_n / n if n > 0 else None)
  row["cys_count"] = str(cys_count)
  row["motif_h3_len_in_8_20"] = "1" if 8 <= n <= 20 else "0"
  row["has_sequence"] = "1"
  return row


def _load_rows(in_csv: Path) -> list[dict[str, str]]:
  if not in_csv.exists():
    return []
  with in_csv.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)


def _build_summary_by_model(rows: list[dict[str, str]]) -> list[dict[str, str]]:
  by_model: dict[str, list[dict[str, str]]] = {}
  for r in rows:
    model = (r.get("model") or "").strip() or "unknown"
    by_model.setdefault(model, []).append(r)

  out: list[dict[str, str]] = []
  for model in sorted(by_model):
    ms = by_model[model]
    n_total = len(ms)
    n_success = sum(1 for r in ms if (r.get("status") or "").strip() == "ok")
    has_structure_n = sum(1 for r in ms if (r.get("has_structure") or "").strip() == "1")
    has_sequence_n = sum(1 for r in ms if (r.get("has_sequence") or "").strip() == "1")
    runtime_values = [_to_float(r.get("duration_sec", "")) for r in ms]
    runtime_values = [x for x in runtime_values if x is not None]
    seq_len_values = [_to_int(r.get("pred_seq_len", "")) for r in ms]
    seq_len_values = [float(x) for x in seq_len_values if x is not None]
    rmsd_values = [_to_float(r.get("cdr_h3_rmsd", "")) for r in ms]
    rmsd_values = [x for x in rmsd_values if x is not None]

    row = {k: "" for k in SUMMARY_BY_MODEL_FIELDS}
    row["model"] = model
    row["n_total"] = str(n_total)
    row["n_success"] = str(n_success)
    row["success_rate"] = _round_or_blank((n_success / n_total) if n_total else None)
    row["has_structure_rate"] = _round_or_blank((has_structure_n / n_total) if n_total else None)
    row["has_sequence_rate"] = _round_or_blank((has_sequence_n / n_total) if n_total else None)
    row["median_runtime_sec"] = _round_or_blank(statistics.median(runtime_values) if runtime_values else None)
    row["median_seq_len"] = _round_or_blank(statistics.median(seq_len_values) if seq_len_values else None)
    row["median_cdr_h3_rmsd"] = _round_or_blank(statistics.median(rmsd_values) if rmsd_values else None)
    row["p25_cdr_h3_rmsd"] = _round_or_blank(_percentile(rmsd_values, 0.25))
    row["p75_cdr_h3_rmsd"] = _round_or_blank(_percentile(rmsd_values, 0.75))
    out.append(row)
  return out


def _build_summary_by_sample(rows: list[dict[str, str]]) -> list[dict[str, str]]:
  by_sample: dict[str, list[dict[str, str]]] = {}
  for r in rows:
    sid = (r.get("sample_id") or "").strip() or "unknown"
    by_sample.setdefault(sid, []).append(r)

  out: list[dict[str, str]] = []
  for sid in sorted(by_sample):
    sr = by_sample[sid]
    n_models_attempted = len({(r.get("model") or "").strip() for r in sr if (r.get("model") or "").strip()})
    n_models_success = len({(r.get("model") or "").strip() for r in sr if (r.get("status") or "").strip() == "ok"})
    best_model = ""
    best_val = None
    for r in sr:
      x = _to_float(r.get("cdr_h3_rmsd", ""))
      if x is None:
        continue
      if best_val is None or x < best_val:
        best_val = x
        best_model = (r.get("model") or "").strip()
    row = {k: "" for k in SUMMARY_BY_SAMPLE_FIELDS}
    row["sample_id"] = sid
    row["n_models_attempted"] = str(n_models_attempted)
    row["n_models_success"] = str(n_models_success)
    row["best_model_by_cdr_h3_rmsd"] = best_model
    row["best_cdr_h3_rmsd"] = _round_or_blank(best_val)
    out.append(row)
  return out


def compute_sequence_metrics(
  in_csv: Path,
  out_csv: Path,
  summary_by_model_csv: Path,
  summary_by_sample_csv: Path,
) -> tuple[Path, Path, Path]:
  rows = _load_rows(in_csv)
  out_rows: list[dict[str, str]] = []
  for row in rows:
    base = {k: (row.get(k, "") or "") for k in EVALUATION_LONG_FIELDS}
    base = _compute_seq_fields(base)
    out_rows.append(base)

  _write_csv(out_csv, EVALUATION_LONG_FIELDS, out_rows)
  _write_csv(summary_by_model_csv, SUMMARY_BY_MODEL_FIELDS, _build_summary_by_model(out_rows))
  _write_csv(summary_by_sample_csv, SUMMARY_BY_SAMPLE_FIELDS, _build_summary_by_sample(out_rows))
  return out_csv, summary_by_model_csv, summary_by_sample_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="计算 Level-0/1 指标并输出模型/样本汇总")
  parser.add_argument(
    "--in-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/evaluation_long.csv"),
    help="输入长表路径",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/evaluation_long_seq.csv"),
    help="输出带序列指标的长表路径",
  )
  parser.add_argument(
    "--summary-by-model",
    type=Path,
    default=Path("outputs/evaluation/all_models/summary_by_model.csv"),
    help="模型汇总输出路径",
  )
  parser.add_argument(
    "--summary-by-sample",
    type=Path,
    default=Path("outputs/evaluation/all_models/summary_by_sample.csv"),
    help="样本汇总输出路径",
  )
  args = parser.parse_args()

  out_csv, s_model, s_sample = compute_sequence_metrics(
    in_csv=args.in_csv,
    out_csv=args.out_csv,
    summary_by_model_csv=args.summary_by_model,
    summary_by_sample_csv=args.summary_by_sample,
  )
  print(f"[OK] 序列指标长表: {out_csv}")
  print(f"[OK] 模型汇总: {s_model}")
  print(f"[OK] 样本汇总: {s_sample}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

