"""Thin CLI for reusable total ion energy inventory analysis."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))


if __name__ == "__main__":
    from energy_analysis.__main__ import main

    raise SystemExit(main())
