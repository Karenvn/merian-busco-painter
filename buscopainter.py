#!/usr/bin/env python3
"""Compatibility wrapper for ``merian-busco-painter paint``."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from merian_busco_painter.cli import paint_main
from merian_busco_painter.painter import *  # noqa: F403


if __name__ == "__main__":
    paint_main()
