#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import socket
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
FAMILY_NAMES = ("dir", "neu")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Hermes MMS cases that have a valid grid file."
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default="diamond",
        help="Case directory or case name under examples/ib/mms. Defaults to diamond.",
    )
    parser.add_argument(
        "--d",
        "--case",
        "--case-dir",
        dest="directory_option",
        default=None,
        help="Case directory or case name under examples/ib/mms. Overrides the positional directory.",
    )
    parser.add_argument(
        "-n",
        "--nproc",
        type=int,
        default=None,
        help=(
            "Number of MPI processes. If omitted, choose a grid-safe per-run "
            "value capped by the detected Slurm allocation or local CPU count."
        ),
    )
    parser.add_argument(
        "--launcher",
        choices=("auto", "mpirun", "srun"),
        default="auto",
        help="MPI launcher. auto uses srun on Perlmutter and mpirun otherwise.",
    )
    parser.add_argument(
        "--auto-nproc",
        action="store_true",
        help=(
            "Choose or cap MPI processes from each gridfile using "
            "floor((nx - 2*MXG)/(2*MXG)), adjusted to divide the physical x cells."
        ),
    )
    parser.add_argument(
        "--max-nproc",
        type=int,
        default=None,
        help=(
            "Maximum MPI processes for automatic sizing. Defaults to the detected "
            "Slurm allocation, or local CPU count outside Slurm."
        ),
    )
    family_group = parser.add_mutually_exclusive_group()
    family_group.add_argument(
        "--dir",
        dest="family",
        action="store_const",
        const="dir",
        help="Run only Dirichlet MMS cases.",
    )
    family_group.add_argument(
        "--neu",
        dest="family",
        action="store_const",
        const="neu",
        help="Run only Neumann MMS cases.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the commands that would be run without launching Hermes.",
    )
    parser.add_argument(
        "--no-prepare",
        action="store_true",
        help="Do not generate/update resolution run directories from dir/neu templates.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Generate/update resolution run directories, then exit without running Hermes.",
    )
    return parser.parse_args()


def resolve_case_dir(directory: str) -> Path:
    candidate = Path(directory).expanduser()
    if candidate.exists():
        return candidate.resolve()

    candidate = HERE / directory
    if candidate.exists():
        return candidate.resolve()

    raise FileNotFoundError(f"Could not find MMS case directory: {directory}")


def git_root() -> Path:
    try:
        output = subprocess.check_output(
            ["git", "-C", str(HERE), "rev-parse", "--show-toplevel"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"Could not find git root from {HERE}") from error
    return Path(output.strip()).resolve()


def parse_bout_inp(path: Path) -> dict[str, dict[str, str]]:
    sections: dict[str, dict[str, str]] = {}
    current: str | None = None
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1].strip()
            sections.setdefault(current, {})
            continue
        if "=" in line and current is not None:
            key, value = [item.strip() for item in line.split("=", 1)]
            sections[current][key] = value.strip().strip('"').strip("'")
    return sections


def mesh_path(run_dir: Path) -> Path | None:
    bout_path = run_dir / "BOUT.inp"
    try:
        mesh_file = parse_bout_inp(bout_path)["mesh"]["file"]
    except KeyError:
        return None

    path = Path(mesh_file).expanduser()
    if path.is_absolute():
        return path
    return (run_dir / path).resolve()


def grid_filename_regex(filename: str) -> re.Pattern[str]:
    match = re.fullmatch(r"(?P<prefix>.*\D)?(?P<size>\d+)(?P<suffix>\D*)", filename)
    if match is None:
        raise ValueError(f"Could not infer a resolution number from grid name: {filename}")
    prefix = match.group("prefix") or ""
    suffix = match.group("suffix") or ""
    return re.compile(rf"^{re.escape(prefix)}(?P<size>\d+){re.escape(suffix)}$")


def template_family_dirs(base: Path, family: str | None) -> list[Path]:
    if base.name in FAMILY_NAMES and (base / "BOUT.inp").is_file():
        return [base] if matches_family(base, family) else []

    dirs: list[Path] = []
    for name in FAMILY_NAMES:
        parent = base / name
        if (parent / "BOUT.inp").is_file() and matches_family(parent, family):
            dirs.append(parent)
    return dirs


def template_gridfiles(template_path: Path) -> list[tuple[int, Path]]:
    mesh_file = parse_bout_inp(template_path)["mesh"]["file"]
    mesh_path = Path(mesh_file).expanduser()
    resolved = mesh_path if mesh_path.is_absolute() else (template_path.parent / mesh_path).resolve()
    grid_dir = resolved.parent
    pattern = grid_filename_regex(resolved.name)
    if not grid_dir.is_dir():
        return []

    grids: list[tuple[int, Path]] = []
    for grid_path in grid_dir.iterdir():
        match = pattern.fullmatch(grid_path.name)
        if match is not None and grid_path.is_file():
            grids.append((int(match.group("size")), grid_path.resolve()))
    return sorted(grids)


def render_template(template_path: Path, mesh_file: str) -> str:
    lines = template_path.read_text().splitlines(keepends=True)
    in_mesh = False
    replaced = False

    for index, line in enumerate(lines):
        stripped = line.split("#", 1)[0].strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_mesh = stripped.lower() == "[mesh]"
            continue
        if in_mesh and re.match(r"\s*file\s*=", line):
            newline = "\r\n" if line.endswith("\r\n") else "\n" if line.endswith("\n") else ""
            lines[index] = f'file = "{mesh_file}"{newline}'
            replaced = True
            break

    if not replaced:
        raise ValueError(f"Could not find [mesh] file entry in {template_path}")
    return "".join(lines)


def prepare_run_dirs(base: Path, family: str | None, dry_run: bool) -> None:
    for parent in template_family_dirs(base, family):
        template_path = parent / "BOUT.inp"
        try:
            grids = template_gridfiles(template_path)
        except (KeyError, ValueError) as error:
            print(f"Skipping template {template_path}: {error}", file=sys.stderr)
            continue

        if not grids:
            print(f"No matching grid files found for template {template_path}", file=sys.stderr)
            continue

        for size, grid_path in grids:
            run_dir = parent / str(size)
            mesh_file = os.path.relpath(grid_path, run_dir)
            target = run_dir / "BOUT.inp"
            text = render_template(template_path, mesh_file)

            if dry_run:
                print(f"Would prepare {target} with mesh {mesh_file}")
                continue

            run_dir.mkdir(parents=True, exist_ok=True)
            if target.exists() and target.read_text() == text:
                continue
            target.write_text(text)


def has_child_cases(path: Path) -> bool:
    return any(child.is_dir() and (child / "BOUT.inp").is_file() for child in path.iterdir())


def natural_key(path: Path) -> list[object]:
    parts = str(path).split("/")
    key: list[object] = []
    for part in parts:
        key.extend(int(piece) if piece.isdigit() else piece for piece in re.split(r"(\d+)", part))
    return key


def candidate_dirs(base: Path) -> list[Path]:
    if (
        (base / "BOUT.inp").is_file()
        and not has_child_cases(base)
        and base.name not in FAMILY_NAMES
    ):
        return [base]

    dirs = []
    for bout_path in base.rglob("BOUT.inp"):
        run_dir = bout_path.parent
        if run_dir.name in FAMILY_NAMES:
            continue
        if not has_child_cases(run_dir):
            dirs.append(run_dir)
    return sorted(dirs, key=natural_key)


def matches_family(run_dir: Path, family: str | None) -> bool:
    if family is None:
        return True
    try:
        parts = run_dir.relative_to(HERE).parts
    except ValueError:
        parts = run_dir.parts
    return family in parts


def runnable_dirs(
    base: Path, family: str | None
) -> tuple[list[tuple[Path, Path]], list[tuple[Path, str]]]:
    runnable: list[tuple[Path, Path]] = []
    skipped: list[tuple[Path, str]] = []

    for run_dir in candidate_dirs(base):
        if not matches_family(run_dir, family):
            continue
        grid_path = mesh_path(run_dir)
        if grid_path is None:
            skipped.append((run_dir, "no [mesh] file in BOUT.inp"))
            continue
        if not grid_path.is_file():
            skipped.append((run_dir, f"missing grid file {grid_path}"))
            continue
        runnable.append((run_dir, grid_path))

    return runnable, skipped


def format_cmd(command: list[str]) -> str:
    return " ".join(command)


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


def selected_launcher(launcher: str) -> str:
    if launcher != "auto":
        return launcher
    return "srun" if is_perlmutter() else "mpirun"


def positive_int(value: str | None) -> int | None:
    if value is None:
        return None
    try:
        parsed = int(value)
    except ValueError:
        return None
    return parsed if parsed > 0 else None


def parse_slurm_cpus_per_node(value: str | None) -> int | None:
    if not value:
        return None

    total = 0
    for item in value.split(","):
        match = re.fullmatch(r"\s*(?P<count>\d+)(?:\(x(?P<repeat>\d+)\))?\s*", item)
        if match is None:
            return None
        count = int(match.group("count"))
        repeat = int(match.group("repeat") or 1)
        total += count * repeat
    return total if total > 0 else None


def detected_max_nproc(override: int | None) -> tuple[int, str]:
    if override is not None:
        return override, "--max-nproc"

    for name in ("SLURM_NTASKS", "SLURM_NPROCS"):
        value = positive_int(os.environ.get(name))
        if value is not None:
            return value, name

    cpus_per_task = positive_int(os.environ.get("SLURM_CPUS_PER_TASK")) or 1
    job_cpus = parse_slurm_cpus_per_node(os.environ.get("SLURM_JOB_CPUS_PER_NODE"))
    if job_cpus is not None:
        return max(1, job_cpus // cpus_per_task), "SLURM_JOB_CPUS_PER_NODE"

    cpus_on_node = positive_int(os.environ.get("SLURM_CPUS_ON_NODE"))
    nodes = positive_int(os.environ.get("SLURM_JOB_NUM_NODES")) or 1
    if cpus_on_node is not None:
        return max(1, (cpus_on_node * nodes) // cpus_per_task), "SLURM_CPUS_ON_NODE"

    return max(1, os.cpu_count() or 1), "os.cpu_count()"


def scalar_grid_int(value: object, name: str, grid_path: Path) -> int:
    import numpy as np

    if value is None:
        raise KeyError(f"{name} was not found in {grid_path}")
    values = np.asarray(value).reshape(-1)
    if values.size == 0:
        raise ValueError(f"{name} in {grid_path} is empty")
    return int(values[0])


def read_grid_decomposition(grid_path: Path) -> tuple[int, int]:
    from boututils.datafile import DataFile

    with DataFile(str(grid_path)) as datafile:
        nx = scalar_grid_int(datafile.read("nx"), "nx", grid_path)
        mxg = scalar_grid_int(datafile.read("MXG"), "MXG", grid_path)
    return nx, mxg


def grid_safe_nproc(grid_path: Path, rank_cap: int) -> tuple[int, str]:
    nx, mxg = read_grid_decomposition(grid_path)
    physical_x = nx - 2 * mxg
    if physical_x <= 0:
        return 1, f"nx={nx}, MXG={mxg}, physical x={physical_x}"

    min_x_per_rank = max(1, 2 * mxg)
    limit = max(1, min(rank_cap, physical_x // min_x_per_rank))
    for candidate in range(limit, 0, -1):
        if physical_x % candidate == 0:
            return (
                candidate,
                f"nx={nx}, MXG={mxg}, physical x={physical_x}, min x/rank={min_x_per_rank}",
            )

    return 1, f"nx={nx}, MXG={mxg}, physical x={physical_x}, min x/rank={min_x_per_rank}"


def choose_nproc(
    args: argparse.Namespace,
    launcher: str,
    grid_path: Path,
    max_nproc: int,
    max_source: str,
) -> tuple[int, str | None]:
    auto_default = args.nproc is None
    if args.auto_nproc or auto_default:
        rank_cap = max_nproc if args.nproc is None else min(args.nproc, max_nproc)
        safe_nproc, reason = grid_safe_nproc(grid_path, rank_cap)
        reason = f"{reason}, rank cap={rank_cap}, system max={max_nproc} from {max_source}"
        if args.nproc is None:
            return safe_nproc, f"auto nproc: {safe_nproc} ({reason})"

        if safe_nproc != args.nproc:
            return safe_nproc, f"auto nproc cap: requested {args.nproc}, using {safe_nproc} ({reason})"
        return safe_nproc, f"auto nproc: using requested {safe_nproc} ({reason})"

    return args.nproc or 1, None


def launcher_command(launcher: str, nproc: int, executable: Path, run_dir: Path) -> list[str]:
    if launcher == "srun":
        return [
            "srun",
            "-n",
            str(nproc),
            str(executable),
            "-d",
            str(run_dir),
        ]

    return [
        "mpirun",
        "-n",
        str(nproc),
        "--use-hwthread-cpus",
        str(executable),
        "-d",
        str(run_dir),
    ]


def main() -> int:
    args = parse_args()
    if args.nproc is not None and args.nproc < 1:
        raise ValueError("--nproc must be at least 1")
    if args.max_nproc is not None and args.max_nproc < 1:
        raise ValueError("--max-nproc must be at least 1")

    base = resolve_case_dir(args.directory_option or args.directory)
    executable = git_root() / "build" / "hermes-3"
    launcher = selected_launcher(args.launcher)
    max_nproc, max_source = detected_max_nproc(args.max_nproc)

    if not args.no_prepare:
        prepare_run_dirs(base, args.family, args.dry_run)
        if args.prepare_only:
            return 0

    runs, skipped = runnable_dirs(base, args.family)

    if not runs:
        family = f" {args.family}" if args.family else ""
        print(f"No runnable{family} cases found under {base}", file=sys.stderr)
        for run_dir, reason in skipped:
            print(f"skipped {run_dir}: {reason}", file=sys.stderr)
        return 1

    if not executable.is_file() and not args.dry_run:
        raise FileNotFoundError(f"Could not find Hermes executable: {executable}")

    for run_dir, grid_path in runs:
        nproc, nproc_note = choose_nproc(args, launcher, grid_path, max_nproc, max_source)
        command = launcher_command(launcher, nproc, executable, run_dir)
        print(f"\nRunning {run_dir}")
        print(f"grid: {grid_path}")
        print(f"launcher: {launcher}")
        if nproc_note:
            print(nproc_note)
        print(format_cmd(command))
        if args.dry_run:
            continue

        completed = subprocess.run(command)
        if completed.returncode != 0:
            return completed.returncode

    if skipped:
        print("\nSkipped cases:")
        for run_dir, reason in skipped:
            print(f"  {run_dir}: {reason}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
