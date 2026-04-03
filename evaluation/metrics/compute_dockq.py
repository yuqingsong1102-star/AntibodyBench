#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path

from evaluation.schema import EVALUATION_LONG_FIELDS


_DOCKQ_RE = re.compile(r"DockQ(?:\s+score)?\s*[:=]\s*([0-9]*\.?[0-9]+)", re.IGNORECASE)
_IRMSD_RE = re.compile(r"iRMSD\s*[:=]\s*([0-9]*\.?[0-9]+)", re.IGNORECASE)
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DOCKQ_CMD = str(ROOT / "scripts" / "tools" / "bin" / "DockQ.py")


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


def _run_dockq(model_path: str, native_path: str, cmd: str = DEFAULT_DOCKQ_CMD) -> tuple[str, str, str]:
  cmd_candidates = [cmd]

  proc = None
  used_cmd = ""
  last_missing = False
  for c in cmd_candidates:
    try:
      proc = subprocess.run(
        [c, model_path, native_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
      )
      used_cmd = c
      break
    except FileNotFoundError:
      last_missing = True
      continue
    except Exception as e:
      return "", "", f"dockq_exec_exception:{type(e).__name__}"

  if proc is None:
    if last_missing:
      return "", "", "dockq_tool_missing"
    return "", "", "dockq_exec_failed"

  text = (proc.stdout or "") + "\n" + (proc.stderr or "")
  m_dockq = _DOCKQ_RE.search(text)
  m_irmsd = _IRMSD_RE.search(text)
  if proc.returncode != 0 and not m_dockq:
    return "", "", f"dockq_failed_rc={proc.returncode}:{used_cmd}"

  dockq_score = ""
  irmsd = ""
  if m_dockq:
    try:
      dockq_score = str(round(float(m_dockq.group(1)), 6))
    except Exception:
      dockq_score = ""
  if m_irmsd:
    try:
      irmsd = str(round(float(m_irmsd.group(1)), 6))
    except Exception:
      irmsd = ""

  if not dockq_score and not irmsd:
    return "", "", "dockq_parse_failed"
  return dockq_score, irmsd, ""


def compute_dockq(in_csv: Path, out_csv: Path, dockq_cmd: str = DEFAULT_DOCKQ_CMD) -> Path:
  rows = _load_rows(in_csv)
  out_rows: list[dict[str, str]] = []

  for row in rows:
    cur = {k: (row.get(k, "") or "") for k in EVALUATION_LONG_FIELDS}
    cur["dockq_score"] = cur.get("dockq_score", "")
    cur["irmsd"] = cur.get("irmsd", "")

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

    dockq_score, irmsd, err = _run_dockq(model_path=model_path, native_path=native_path, cmd=dockq_cmd)
    if dockq_score:
      cur["dockq_score"] = dockq_score
    if irmsd:
      cur["irmsd"] = irmsd
    if err:
      cur["metric_error_code"] = err
      cur["metric_error"] = err
    out_rows.append(cur)

  _write_rows(out_csv, out_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="批量计算 DockQ / iRMSD 并回写长表")
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
  parser.add_argument("--dockq-cmd", default=DEFAULT_DOCKQ_CMD, help="DockQ 命令路径")
  args = parser.parse_args()

  out = compute_dockq(in_csv=args.in_csv, out_csv=args.out_csv, dockq_cmd=args.dockq_cmd)
  print(f"[OK] 写出 DockQ 指标: {out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

