#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path

from evaluation.schema import EVALUATION_LONG_FIELDS


_TM_RE = re.compile(r"TM-score\s*=\s*([0-9]*\.?[0-9]+)")
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TMSCORE_CMD = str(ROOT / "scripts" / "tools" / "bin" / "TMscore")


def _load_rows(in_csv: Path) -> list[dict[str, str]]:
  if not in_csv.exists():
    return []
  with in_csv.open("r", encoding="utf-8", newline="") as f:
    return list(csv.DictReader(f))


def _write_rows(out_csv: Path, rows: list[dict[str, str]]) -> None:
  out_csv.parent.mkdir(parents=True, exist_ok=True)
  with out_csv.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=EVALUATION_LONG_FIELDS)
    writer.writeheader()
    writer.writerows(rows)


def _run_tmscore(model_path: str, native_path: str, cmd: str = DEFAULT_TMSCORE_CMD) -> tuple[str, str]:
  try:
    proc = subprocess.run(
      [cmd, model_path, native_path],
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      text=True,
      check=False,
    )
  except FileNotFoundError:
    return "", "tmscore_tool_missing"
  except Exception as e:
    return "", f"tmscore_exec_exception:{type(e).__name__}"

  text = (proc.stdout or "") + "\n" + (proc.stderr or "")
  m = _TM_RE.search(text)
  if not m:
    if proc.returncode != 0:
      return "", f"tmscore_failed_rc={proc.returncode}"
    return "", "tmscore_parse_failed"
  try:
    return str(round(float(m.group(1)), 6)), ""
  except Exception:
    return "", "tmscore_parse_failed"


def compute_tm_score(in_csv: Path, out_csv: Path, tmscore_cmd: str = DEFAULT_TMSCORE_CMD) -> Path:
  rows = _load_rows(in_csv)
  out_rows: list[dict[str, str]] = []

  for row in rows:
    cur = {k: (row.get(k, "") or "") for k in EVALUATION_LONG_FIELDS}
    cur["tm_score"] = cur.get("tm_score", "")

    model_path = (cur.get("structure_path") or "").strip()
    native_path = (cur.get("reference_complex_path") or "").strip()

    if not model_path:
      cur["metric_error_code"] = "missing_pred_structure"
      out_rows.append(cur)
      continue
    if not native_path:
      cur["metric_error_code"] = "missing_reference_complex"
      out_rows.append(cur)
      continue

    score, err = _run_tmscore(model_path=model_path, native_path=native_path, cmd=tmscore_cmd)
    if score:
      cur["tm_score"] = score
    if err:
      cur["metric_error_code"] = err
      cur["metric_error"] = err
    out_rows.append(cur)

  _write_rows(out_csv, out_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="批量计算 TM-score 并回写长表")
  parser.add_argument(
    "--in-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/evaluation_long_metrics.csv"),
    help="输入长表",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/all_models/evaluation_long_metrics.csv"),
    help="输出长表",
  )
  parser.add_argument("--tmscore-cmd", default=DEFAULT_TMSCORE_CMD, help="TMscore 命令路径")
  args = parser.parse_args()

  out = compute_tm_score(in_csv=args.in_csv, out_csv=args.out_csv, tmscore_cmd=args.tmscore_cmd)
  print(f"[OK] 写出 TM-score 指标: {out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

