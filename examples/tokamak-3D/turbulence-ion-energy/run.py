#!/usr/bin/env python3
"""Launch this example into a fresh output directory."""

import argparse
import shlex
import shutil
from pathlib import Path

from boututils.run_wrapper import launch_safe

here = Path(__file__).resolve().parent
cli = argparse.ArgumentParser(description=__doc__)
cli.add_argument("directory", type=Path)
cli.add_argument("--input", type=Path, default=here / "BOUT.inp")
cli.add_argument("--exe", type=Path, default=here.parents[2] / "build/hermes-3")
cli.add_argument(
    "--grid",
    type=Path,
    default=here.parents[2] / "build/examples/tokamak-3D/tokamak.nc",
)
cli.add_argument("--nproc", type=int, default=6)
cli.add_argument("--mz", type=int, default=1)
cli.add_argument("--nout", type=int, default=10)
cli.add_argument("--timestep", type=float, default=1.0)
cli.add_argument("--diagnose", action="store_true")
cli.add_argument("--set", action="append", default=[], metavar="OPTION=VALUE")
cli.add_argument("--restart-from", type=Path)
args = cli.parse_args()
args.directory.mkdir(parents=True, exist_ok=False)
shutil.copyfile(args.input, args.directory / "BOUT.inp")
if args.restart_from:
    for path in args.restart_from.glob("BOUT.restart.*"):
        shutil.copyfile(path, args.directory / path.name)
command = [
    str(args.exe.resolve()),
    "-d",
    str(args.directory.resolve()),
    "mesh:file=" + str(args.grid.resolve()),
    f"MZ={args.mz}",
    f"nout={args.nout}",
    f"timestep={args.timestep}",
]
if args.diagnose:
    command += [
        f"{name}:diagnose=true" for name in ("e", "d+", "evolve_vorticity_ion_energy")
    ]
if args.restart_from:
    command += ["restart=true", "append=false"]
command += args.set
launch_safe(
    shlex.join(command),
    nproc=args.nproc,
    mthread=1,
    output=shlex.quote(str(args.directory / "run.log")),
)
