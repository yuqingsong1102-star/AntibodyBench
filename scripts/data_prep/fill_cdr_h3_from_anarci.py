#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from runpy import run_path


def main() -> int:
  script_dir = Path(__file__).resolve().parent
  target = script_dir.parent / "fill_cdr_h3_from_anarci.py"
  run_path(str(target), run_name="__main__")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

