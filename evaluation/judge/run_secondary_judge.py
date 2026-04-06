#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
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


def run_secondary_judge(seed_level_csv: Path, out_csv: Path, judge_manifest: Path | None = None) -> Path:
  rows = _load_csv(seed_level_csv)
  secondary_by_candidate_seed: dict[tuple[str, str], dict[str, str]] = {}
  if judge_manifest and judge_manifest.exists():
    for row in _load_csv(judge_manifest):
      key = ((row.get("candidate_id") or "").strip(), (row.get("judge_seed") or "0").strip())
      if key[0]:
        secondary_by_candidate_seed[key] = row

  out_rows: list[dict[str, str]] = []
  for row in rows:
    cur = {k: (row.get(k, "") or "") for k in EXTERNAL_SEED_FIELDS}
    key = ((cur.get("candidate_id") or "").strip(), (cur.get("judge_seed") or "0").strip())
    imported = secondary_by_candidate_seed.get(key)
    if imported:
      # 仅记录 agreement 标志，不覆盖 primary judge 主表。
      primary_ok = (cur.get("external_status") or "") == "ok"
      secondary_ok = ((imported.get("external_status") or imported.get("status") or "") == "ok")
      cur["secondary_judge_agreement"] = "1" if primary_ok == secondary_ok else "0"
    elif not cur.get("secondary_judge_agreement"):
      cur["secondary_judge_agreement"] = ""
    out_rows.append(cur)

  _write_csv(out_csv, out_rows)
  return out_csv


def main() -> int:
  parser = argparse.ArgumentParser(description="导入 secondary judge 一致性结果（可选）")
  parser.add_argument(
    "--seed-level-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="seed 级主表",
  )
  parser.add_argument(
    "--out-csv",
    type=Path,
    default=Path("outputs/evaluation/external_binding_benchmark/external_eval_seed_level.csv"),
    help="输出 seed 级主表（可与输入相同）",
  )
  parser.add_argument("--judge-manifest", type=Path, default=None, help="secondary judge 结果 manifest")
  args = parser.parse_args()
  out = run_secondary_judge(args.seed_level_csv, args.out_csv, judge_manifest=args.judge_manifest)
  print(f"[OK] secondary judge 合并: {out}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
