#!/usr/bin/env python3
"""Run a total-ion-energy resolution scan; analyse_scan.py reads the outputs."""

import argparse
import shlex
import shutil
from pathlib import Path

from boututils.run_wrapper import launch_safe

HERE = Path(__file__).resolve().parent
cli = argparse.ArgumentParser(description=__doc__)
cli.add_argument("directory", type=Path)
cli.add_argument("--exe", type=Path, default=HERE.parents[2] / "build/hermes-3")
cli.add_argument("--levels", type=int, default=3)
cli.add_argument("--base-nx", type=int, default=32)
cli.add_argument("--base-nz", type=int, default=32)
cli.add_argument("--nproc", type=int, default=1)
cli.add_argument("--end-time", type=float, default=20000)
cli.add_argument("--outputs", type=int, default=200)
cli.add_argument("--set", action="append", default=[], metavar="OPTION=VALUE")
args = cli.parse_args()
if min(args.levels, args.base_nx, args.base_nz, args.nproc, args.outputs) < 1:
    cli.error("Grid sizes, ranks, levels and outputs must be positive")
args.directory.mkdir(parents=True, exist_ok=False)
for level in range(args.levels):
    nx, nz = args.base_nx * 2**level, args.base_nz * 2**level
    output = args.directory / f"{nx}x{nz}"
    output.mkdir()
    shutil.copyfile(HERE / "BOUT.inp", output / "BOUT.inp")
    command = [
        str(args.exe.resolve()),
        "-d",
        str(output.resolve()),
        f"mesh:nx={nx + 4}",
        f"mesh:nz={nz}",
        f"NXPE={args.nproc}",
        f"solver:nout={args.outputs}",
        f"solver:output_step={args.end_time / args.outputs}",
        *args.set,
    ]
    launch_safe(
        shlex.join(command),
        nproc=args.nproc,
        mthread=1,
        output=shlex.quote(str(output / "run.log")),
    )
