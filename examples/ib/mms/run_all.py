#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shlex
import socket
import subprocess
import sys
from pathlib import Path

import grid_gen


HERE = Path(__file__).resolve().parent
DEFAULT_GRID_NP = 64
DEFAULT_GRID_LEVELS = 4
PERLMUTTER_GRID_LEVELS = 6
LOCAL_DEFAULT_SKIP_CASES = {"d3d"}
CASE_DEFAULTS = {
    "d3d": {
        "np": 256,
        "levels": 4,
        "interp": True,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate grids, run Hermes-3, and plot MMS diagnostics."
    )
    parser.add_argument(
        "--d",
        "--case",
        dest="cases",
        action="append",
        default=None,
        help=(
            "MMS case to process, such as circle or diamond. Repeat for multiple "
            "cases. Defaults to all cases supported by grid_gen.py."
        ),
    )
    parser.add_argument(
        "--np",
        type=int,
        default=None,
        metavar="NPOINTS",
        help="Starting perpendicular gridpoints for grid_gen.py. Defaults to 64, or 256 for d3d.",
    )
    parser.add_argument(
        "--levels",
        "--nlevels",
        "--grid-n",
        dest="levels",
        type=int,
        default=None,
        metavar="NFILES",
        help=(
            "Number of resolution levels for grid_gen.py. Defaults to 4, "
            "or 6 on Perlmutter for regular cases. d3d defaults to 4."
        ),
    )
    parser.add_argument(
        "-n",
        "--nproc",
        type=int,
        default=None,
        help="Number of MPI processes for run.py. Defaults to run.py's host-aware choice.",
    )
    parser.add_argument(
        "--launcher",
        choices=("auto", "mpirun", "srun"),
        default="auto",
        help="MPI launcher passed to run.py. auto uses srun on Perlmutter and mpirun otherwise.",
    )
    parser.add_argument(
        "--auto-nproc",
        action="store_true",
        help="Pass --auto-nproc to run.py to choose/cap ranks from each gridfile.",
    )
    parser.add_argument(
        "--max-nproc",
        type=int,
        default=None,
        help="Maximum MPI processes for run.py automatic sizing.",
    )
    parser.add_argument(
        "--tind",
        type=int,
        default=None,
        help="Time index passed to mms.py. Defaults to mms.py's halfway timestep.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("out"),
        help="Directory for generated PDFs. Defaults to ./out.",
    )
    parser.add_argument(
        "--force-grids",
        action="store_true",
        help="Regenerate existing grid files.",
    )
    parser.add_argument(
        "--skip-grid",
        action="store_true",
        help="Skip grid generation.",
    )
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Skip Hermes runs.",
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Skip MMS plotting.",
    )
    parser.add_argument(
        "--plot",
        "--show",
        dest="plot",
        action="store_true",
        help="Show mms.py plots interactively after saving.",
    )
    parser.add_argument(
        "--interp",
        action="store_true",
        help="Pass --interp to mms.py for every selected case.",
    )
    parser.add_argument(
        "--log",
        dest="log_error",
        action="store_true",
        help="Pass --log to mms.py so spatial error panels use log10(error).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the workflow commands without running them.",
    )
    return parser.parse_args()


def selected_cases(args: argparse.Namespace) -> list[str]:
    if args.cases:
        return args.cases
    cases = set(grid_gen.CASES)
    for child in HERE.iterdir():
        if child.is_dir() and ((child / "dir").is_dir() or (child / "neu").is_dir()):
            cases.add(child.name)
    if not is_perlmutter():
        cases.difference_update(LOCAL_DEFAULT_SKIP_CASES)
    return sorted(cases)


def case_np(case: str, args: argparse.Namespace) -> int:
    if args.np is not None:
        return args.np
    return int(CASE_DEFAULTS.get(case, {}).get("np", DEFAULT_GRID_NP))


def host_values() -> tuple[str, ...]:
    values = []
    for key in ("NERSC_HOST", "HOST", "HOSTNAME"):
        value = os.environ.get(key)
        if value:
            values.append(value)
    values.append(socket.gethostname())
    return tuple(value.lower() for value in values if value)


def is_perlmutter() -> bool:
    return any(value == "perlmutter" or "perlmutter" in value for value in host_values())


def case_levels(case: str, args: argparse.Namespace) -> int:
    if args.levels is not None:
        return args.levels
    if "levels" in CASE_DEFAULTS.get(case, {}):
        return int(CASE_DEFAULTS[case]["levels"])
    if is_perlmutter():
        return PERLMUTTER_GRID_LEVELS
    return DEFAULT_GRID_LEVELS


def case_uses_interp(case: str, args: argparse.Namespace) -> bool:
    return args.interp or bool(CASE_DEFAULTS.get(case, {}).get("interp", False))


def command_string(command: list[str]) -> str:
    return shlex.join(command)


def run_command(command: list[str], dry_run: bool) -> None:
    print(command_string(command))
    if dry_run:
        return
    completed = subprocess.run(command, cwd=HERE)
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)


def grid_command(case: str, args: argparse.Namespace) -> list[str] | None:
    if case not in grid_gen.CASES:
        print(
            f"Skipping grid generation for {case}: not configured in grid_gen.py "
            f"(default would be --np {case_np(case, args)} --n {case_levels(case, args)})"
        )
        return None

    command = [
        sys.executable,
        "-B",
        str(HERE / "grid_gen.py"),
        "--d",
        case,
        "--np",
        str(case_np(case, args)),
        "--n",
        str(case_levels(case, args)),
    ]
    if args.force_grids:
        command.append("--force")
    if args.dry_run:
        command.append("--dry-run")
    return command


def run_case_command(case: str, args: argparse.Namespace) -> list[str]:
    command = [
        sys.executable,
        "-B",
        str(HERE / "run.py"),
        "--d",
        case,
    ]
    if args.nproc is not None:
        command.extend(["-n", str(args.nproc)])
    if args.launcher != "auto":
        command.extend(["--launcher", args.launcher])
    if args.auto_nproc:
        command.append("--auto-nproc")
    if args.max_nproc is not None:
        command.extend(["--max-nproc", str(args.max_nproc)])
    if args.dry_run:
        command.append("--dry-run")
    return command


def has_bout_inputs(case: str) -> bool:
    case_dir = HERE / case
    return case_dir.is_dir() and any(case_dir.rglob("BOUT.inp"))


def has_completed_runs(case: str) -> bool:
    case_dir = HERE / case
    return case_dir.is_dir() and any(case_dir.rglob("BOUT.dmp.*.nc"))


def plot_commands(case: str, args: argparse.Namespace) -> list[list[str]]:
    commands: list[list[str]] = []
    use_interp = case_uses_interp(case, args)
    mode_suffix = "Interp" if use_interp else ""
    for family, extra in (("Neu", []), ("Dir", ["--dir"])):
        command = [
            sys.executable,
            "-B",
            str(HERE / "mms.py"),
            "--d",
            case,
            "--out-dir",
            str(args.out_dir),
            "--output",
            str(args.out_dir / f"{case}{family}{mode_suffix}MMS.pdf"),
            *extra,
        ]
        if use_interp:
            command.append("--interp")
        else:
            command.extend(
                [
                    "--ghost-output",
                    str(args.out_dir / f"{case}{family}GhostMMS.pdf"),
                ]
            )
        if args.tind is not None:
            command.extend(["--tind", str(args.tind)])
        if args.plot:
            command.append("--plot")
        if args.log_error:
            command.append("--log")
        commands.append(command)
    return commands


def main() -> int:
    args = parse_args()
    if args.np is not None and args.np < 1:
        raise SystemExit("--np must be >= 1")
    if args.levels is not None and args.levels < 1:
        raise SystemExit("--levels must be >= 1")
    if args.nproc is not None and args.nproc < 1:
        raise SystemExit("--nproc must be >= 1")
    if args.max_nproc is not None and args.max_nproc < 1:
        raise SystemExit("--max-nproc must be >= 1")

    for case in selected_cases(args):
        case_dir = HERE / case
        if case not in grid_gen.CASES and not case_dir.is_dir():
            raise FileNotFoundError(f"Could not find MMS case directory: {case}")

        print(f"\n== {case} ==")

        if not args.skip_grid:
            command = grid_command(case, args)
            if command is not None:
                run_command(command, args.dry_run)

        if not args.skip_run:
            if args.dry_run or has_bout_inputs(case):
                run_command(run_case_command(case, args), args.dry_run)
            else:
                print(f"Skipping runs for {case}: no BOUT.inp templates or run directories found.")

        if not args.skip_plot:
            if args.dry_run or has_completed_runs(case):
                for command in plot_commands(case, args):
                    run_command(command, args.dry_run)
            else:
                print(f"Skipping plots for {case}: no BOUT.dmp.*.nc files found.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
