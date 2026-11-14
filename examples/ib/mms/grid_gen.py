#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import os
from types import SimpleNamespace
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent

@dataclass(frozen=True)
class GridCase:
    module: str
    prefix: str
    options: dict[str, object]

    def filename(self, size: int) -> str:
        return f"{self.prefix}{size}.fci.nc"


CASES = {
    "diamond": (
        GridCase(
            module="zoidberg.ib.mms_square",
            prefix="sqrt",
            options={"square_rot": True},
        ),
    ),
    "circle": (
        GridCase(
            module="zoidberg.ib.mms_circle",
            prefix="cl",
            options={},
        ),
        GridCase(
            module="zoidberg.ib.mms_circle",
            prefix="clShift",
            options={"x0": -1.0, "x1": 1.0, "z0": -1.0, "z1": 1.0},
        ),
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate IB MMS FCI grid files using the Zoidberg IB scripts."
    )
    parser.add_argument(
        "--d",
        "--case",
        dest="case",
        choices=sorted(CASES),
        required=True,
        help="MMS geometry/case to generate.",
    )

    parser.add_argument(
        "--n",
        type=int,
        default=4,
        metavar="NFILES",
        help="Number of resolution levels to generate.",
    )
    parser.add_argument(
        "--np",
        type=int,
        default=64,
        metavar="NPOINTS",
        help="Starting perpendicular gridpoints in x and z.",
    )

    parser.add_argument("--force", action="store_true", help="Overwrite existing grids.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print Zoidberg commands without running them.",
    )
    return parser.parse_args()


def selected_sizes(args: argparse.Namespace) -> tuple[int, ...]:
    if args.n < 1:
        raise SystemExit("--n must be >= 1")
    if args.np < 1:
        raise SystemExit("--np must be >= 1")
    if args.np & (args.np - 1):
        raise SystemExit("--np must be a power of two")

    return tuple(args.np * 2**i for i in range(args.n))


def setup_plot_environment() -> None:
    os.environ.setdefault("MPLBACKEND", "Agg")
    mplconfigdir = Path(os.environ.get("TMPDIR", "/tmp")) / "hermes_ib_mplconfig"
    mplconfigdir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mplconfigdir))


def zoidberg_args(grid_case: GridCase, size: int, target: Path) -> SimpleNamespace:
    values = {
        "nx": size,
        "nz": size,
        "ny": 1,
        "Ly": 1.0,
        "x0": 0.0,
        "x1": 1.0,
        "z0": 0.0,
        "z1": 1.0,
        "mxg": 4,
        "circle_n": None,
        "square_mode": "coords",
        "square_rot": False,
        "nsamples_per_dim": 1,
        "partial_boundaries": True,
        "fname": str(target),
        "plot": False,
    }
    values.update(grid_case.options)
    return SimpleNamespace(**values)


def describe_call(module: str, args: SimpleNamespace) -> str:
    parts = [
        f"n={args.nx}",
        f"fname={args.fname}",
        "plot=False",
    ]
    if args.square_rot:
        parts.append("rot=True")
    if (args.x0, args.x1, args.z0, args.z1) != (0.0, 1.0, 0.0, 1.0):
        parts.append(f"x=[{args.x0}, {args.x1}]")
        parts.append(f"z=[{args.z0}, {args.z1}]")
    return f"{module}.main(" + ", ".join(parts) + ")"


def generate_one(
    args: argparse.Namespace,
    grid_case: GridCase,
    size: int,
    outdir: Path,
) -> None:
    module = grid_case.module
    target = outdir / grid_case.filename(size)
    if target.exists() and not args.force:
        print(f"skip existing {target}")
        return
    if target.exists():
        target.unlink()

    zargs = zoidberg_args(grid_case, size, target)

    print(describe_call(module, zargs))
    if args.dry_run:
        return

    setup_plot_environment()
    zmodule = importlib.import_module(module)
    zmodule.main(zargs)
    if not target.is_file():
        raise SystemExit(f"Expected Zoidberg output was not created: {target}")


def main() -> None:
    args = parse_args()
    outdir = HERE / args.case / "gridfiles"
    outdir.mkdir(parents=True, exist_ok=True)

    for size in selected_sizes(args):
        for grid_case in CASES[args.case]:
            generate_one(args, grid_case, size, outdir)


if __name__ == "__main__":
    main()
