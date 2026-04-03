#!/usr/bin/env python3
from __future__ import annotations

import runpy
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
TOOL_ROOT = SCRIPT_DIR.parent
VENDOR_DIR = TOOL_ROOT / "vendor"


def main() -> int:
	if not (VENDOR_DIR / "DockQ").exists():
		print(f"[ERROR] Project-local DockQ package not found under {VENDOR_DIR / 'DockQ'}", file=sys.stderr)
		return 127

	sys.path.insert(0, str(VENDOR_DIR))
	sys.argv[0] = str(Path(__file__).resolve())
	runpy.run_module("DockQ", run_name="__main__")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
