"""antibench/__main__.py — 支持 python -m antibench <stage> 用法。"""
import sys
from antibench.cli import main

sys.exit(main())
