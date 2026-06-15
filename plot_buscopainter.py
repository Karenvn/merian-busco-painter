#!/usr/bin/env python3
"""Compatibility wrapper for ``merian-busco-painter plot``."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from merian_busco_painter.cli import plot_main
from merian_busco_painter.plotter import *  # noqa: F403


if __name__ == "__main__":
    plot_main()
