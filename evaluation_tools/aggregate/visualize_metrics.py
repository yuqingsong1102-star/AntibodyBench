#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
  parser = argparse.ArgumentParser(description="简单指标可视化占位脚本")
  parser.add_argument("--model", required=True)
  parser.add_argument(
    "--root",
    type=Path,
    default=Path(__file__).resolve().parents[1],
  )
  args = parser.parse_args()

  ref_csv = args.root / "outputs" / "evaluation" / args.model / "prediction_reference.csv"
  if not ref_csv.exists():
    print(f"[WARN] 文件不存在: {ref_csv}")
    return 0
  with ref_csv.open("r", encoding="utf-8", newline="") as f:
    n = sum(1 for _ in csv.DictReader(f))
  print(f"[OK] {args.model} 当前可视化输入样本数: {n}")
  print("[INFO] 如需图像，请后续接入 matplotlib 绘图逻辑。")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
