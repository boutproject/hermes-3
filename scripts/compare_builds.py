#!/usr/bin/env python3
"""
Compare two hermes-3 builds for output correctness and wall-clock timing.

Runs every integrated test case (or a chosen subset) with both binaries,
then reports per-variable field differences and wall-clock speedup.

Usage (on HPC, after building both binaries):

    python3 scripts/compare_builds.py \\
        --old  /path/to/old/hermes-3 \\
        --new  /path/to/new/hermes-3 \\
        --tests-dir tests/integrated \\
        --mpirun "mpirun -np"

Or for a specific subset of tests:

    python3 scripts/compare_builds.py \\
        --old  /path/to/old/hermes-3 \\
        --new  /path/to/new/hermes-3 \\
        --test tests/integrated/1D-recycling \\
        --test tests/integrated/2D-production \\
        --mpirun "mpirun -np"

Notes
-----
* 2D tests (2D-production, 2D-recycling) download a ~300 MB zip from Zenodo
  the first time they run.  That download is done once by running the original
  runtest; after that, compare_builds.py reuses the cached grid/restart files.
  If those files are missing the test is skipped with a clear warning.
* nproc is auto-detected from each test's runtest script.
* The --mpirun flag should match your cluster's launcher, e.g.:
    Imperial CX3:  "mpirun -np"   (default)
    Slurm systems: "srun -n"
"""

import argparse
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import time

import numpy
from numpy.testing._private.utils import nulp_diff

try:
    import xhermes
except ImportError:
    sys.exit("xhermes is required: pip install xhermes")


# ─── helpers ─────────────────────────────────────────────────────────────────

def detect_nproc(test_dir: pathlib.Path) -> int:
    """Parse the nproc value from the test's runtest script (default 1)."""
    runtest = test_dir / "runtest"
    if not runtest.is_file():
        return 1
    m = re.search(r'nproc\s*=\s*(\d+)', runtest.read_text())
    return int(m.group(1)) if m else 1


def missing_input_files(test_dir: pathlib.Path) -> list[str]:
    """
    Return any .nc files referenced by BOUT.inp that don't exist in test_dir
    or test_dir/data/.  These are typically grid/restart files that the
    runtest downloads from Zenodo.
    """
    bout_inp = test_dir / "data" / "BOUT.inp"
    if not bout_inp.is_file():
        return []
    missing = []
    for m in re.finditer(r'file\s*=\s*"?([^"\s]+\.nc)"?', bout_inp.read_text()):
        fname = m.group(1)
        if not (test_dir / fname).is_file() and not (test_dir / "data" / fname).is_file():
            missing.append(fname)
    # Also check for restart files if the runtest mentions them
    runtest = test_dir / "runtest"
    if runtest.is_file() and "BOUT.restart" in runtest.read_text():
        if not list(test_dir.glob("BOUT.restart.*.nc")):
            missing.append("BOUT.restart.*.nc")
    return missing


def setup_work_dir(test_dir: pathlib.Path,
                   work_dir: pathlib.Path,
                   executable: pathlib.Path) -> None:
    """
    Prepare an isolated working directory for one run:
      work_dir/
        hermes-3          -> symlink to executable
        data/             -> copy of test_dir/data/ (minus any existing dmp files)
        *.nc              -> symlinks to grid/restart .nc files in test_dir
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Symlink the binary
    (work_dir / "hermes-3").symlink_to(executable.resolve())

    # Copy the input directory (BOUT.inp and any other small config files)
    data_dst = work_dir / "data"
    shutil.copytree(test_dir / "data", data_dst)

    # Remove any pre-existing dump files from the copy
    for f in data_dst.glob("BOUT.dmp*.nc"):
        f.unlink()

    # Symlink large read-only inputs (grid files, restart files) from test_dir
    for nc_file in test_dir.glob("*.nc"):
        (work_dir / nc_file.name).symlink_to(nc_file.resolve())


# ─── runner ──────────────────────────────────────────────────────────────────

def run_hermes(work_dir: pathlib.Path,
               nproc: int,
               mpirun: str) -> float:
    """Run hermes-3 from work_dir and return wall-clock seconds."""
    if nproc > 1:
        cmd = mpirun.split() + [str(nproc), "./hermes-3", "-d", "data"]
    else:
        cmd = ["./hermes-3", "-d", "data"]

    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
    elapsed = time.perf_counter() - t0

    if result.returncode != 0:
        print("  [stdout tail]")
        print(result.stdout[-2000:])
        print("  [stderr tail]")
        print(result.stderr[-1000:])
        raise RuntimeError(f"hermes-3 exited with code {result.returncode}")

    return elapsed


# ─── comparison ──────────────────────────────────────────────────────────────

def compare_outputs(data_old: pathlib.Path, data_new: pathlib.Path) -> bool:
    """Compare all field variables in the last timestep. Returns True if ok."""
    try:
        ds_old = xhermes.open(data_old, unnormalise=False)
        ds_new = xhermes.open(data_new, unnormalise=False)
    except Exception as e:
        print(f"  [!] Could not load output: {e}")
        return False

    old_last = ds_old.isel(t=-1)
    new_last = ds_new.isel(t=-1)

    vars_old = set(old_last.data_vars)
    vars_new = set(new_last.data_vars)
    common = sorted(vars_old & vars_new)

    if vars_old - vars_new:
        print(f"  [!] Variables only in OLD: {sorted(vars_old - vars_new)}")
    if vars_new - vars_old:
        print(f"  [!] Variables only in NEW: {sorted(vars_new - vars_old)}")

    print(f"\n  {'Variable':<30} {'max |rel diff|':>16}  {'RMS rel diff':>14}  {'max ULPs':>10}")
    print("  " + "-" * 74)

    any_large_diff = False
    for v in common:
        a = old_last[v].values.ravel().astype(float)
        b = new_last[v].values.ravel().astype(float)

        mask = numpy.isfinite(a) & numpy.isfinite(b) & (a != 0)
        if not mask.any():
            continue

        rel = numpy.abs((b[mask] - a[mask]) / a[mask])
        max_rel = rel.max()
        rms_rel = numpy.sqrt((rel**2).mean())
        max_ulps = int(nulp_diff(a[mask], b[mask]).max())

        flag = "  <-- LARGE" if max_rel > 1e-10 else ""
        if max_rel > 1e-10:
            any_large_diff = True

        print(f"  {v:<30} {max_rel:>16.3e}  {rms_rel:>14.3e}  {max_ulps:>10}{flag}")

    return not any_large_diff


# ─── per-test driver ─────────────────────────────────────────────────────────

def run_test_case(test_dir: pathlib.Path,
                  exec_old: pathlib.Path,
                  exec_new: pathlib.Path,
                  nruns: int,
                  mpirun: str) -> bool:
    name = test_dir.name
    print(f"\n{'='*70}")
    print(f"Test: {name}  (nproc={detect_nproc(test_dir)})")
    print('='*70)

    if not (test_dir / "data" / "BOUT.inp").is_file():
        print(f"  [!] No data/BOUT.inp — skipping")
        return True

    missing = missing_input_files(test_dir)
    if missing:
        print(f"  [!] Skipping: missing input file(s): {missing}")
        print(f"       Run the original runtest once to download them:")
        print(f"       python3 {test_dir}/runtest")
        return True

    nproc = detect_nproc(test_dir)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        work_old = tmp / "old"
        work_new = tmp / "new"

        setup_work_dir(test_dir, work_old, exec_old)
        setup_work_dir(test_dir, work_new, exec_new)

        # ── old binary ───────────────────────────────────────────────────────
        times_old = []
        try:
            for i in range(nruns):
                t = run_hermes(work_old, nproc, mpirun)
                times_old.append(t)
                if i < nruns - 1:
                    for f in (work_old / "data").glob("BOUT.dmp*.nc"):
                        f.unlink()
        except RuntimeError as e:
            print(f"  [!] OLD binary failed: {e}")
            return False

        # ── new binary ───────────────────────────────────────────────────────
        times_new = []
        try:
            for i in range(nruns):
                t = run_hermes(work_new, nproc, mpirun)
                times_new.append(t)
                if i < nruns - 1:
                    for f in (work_new / "data").glob("BOUT.dmp*.nc"):
                        f.unlink()
        except RuntimeError as e:
            print(f"  [!] NEW binary failed: {e}")
            return False

        # ── timing ───────────────────────────────────────────────────────────
        mean_old = numpy.mean(times_old)
        mean_new = numpy.mean(times_new)
        speedup = mean_old / mean_new
        sign = "faster" if speedup > 1 else "slower"

        print(f"\n  Timing ({nruns} run(s) each):")
        print(f"    OLD: {mean_old:.2f} s   ({', '.join(f'{t:.2f}' for t in times_old)})")
        print(f"    NEW: {mean_new:.2f} s   ({', '.join(f'{t:.2f}' for t in times_new)})")
        print(f"    Speedup: {speedup:.3f}x ({abs(speedup-1)*100:.1f}% {sign})")

        # ── correctness ──────────────────────────────────────────────────────
        print("\n  Field differences (last timestep):")
        ok = compare_outputs(work_old / "data", work_new / "data")
        if ok:
            print("\n  => All variables match within 1e-10 relative tolerance.")
        else:
            print("\n  => WARNING: large differences detected (see above).")

        return ok


# ─── main ────────────────────────────────────────────────────────────────────

def discover_tests(tests_dir: pathlib.Path) -> list[pathlib.Path]:
    """Return all subdirectories of tests_dir that contain a runtest script."""
    return sorted(p.parent for p in tests_dir.glob("*/runtest"))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--old", required=True, type=pathlib.Path,
                        help="OLD hermes-3 executable (baseline)")
    parser.add_argument("--new", required=True, type=pathlib.Path,
                        help="NEW hermes-3 executable (with changes)")

    test_group = parser.add_mutually_exclusive_group(required=True)
    test_group.add_argument("--test", action="append", dest="tests",
                            type=pathlib.Path, metavar="DIR",
                            help="Integrated test directory (repeatable)")
    test_group.add_argument("--tests-dir", type=pathlib.Path, metavar="DIR",
                            help="Run all tests found under this directory "
                                 "(e.g. tests/integrated)")

    parser.add_argument("--nruns", type=int, default=3,
                        help="Timed runs per binary per test (default: 3)")
    parser.add_argument("--mpirun", default="mpirun -np",
                        help='MPI launch prefix including the -n flag '
                             '(default: "mpirun -np"). '
                             'Use "srun -n" on Slurm clusters.')
    args = parser.parse_args()

    for exe in (args.old, args.new):
        if not exe.is_file():
            sys.exit(f"Executable not found: {exe}")

    if args.tests_dir:
        tests = discover_tests(args.tests_dir.resolve())
        if not tests:
            sys.exit(f"No runtest scripts found under {args.tests_dir}")
        print(f"Found {len(tests)} test(s) under {args.tests_dir}")
    else:
        tests = [t.resolve() for t in args.tests]

    passed = []
    failed = []
    skipped = []

    for test in tests:
        try:
            ok = run_test_case(test, args.old.resolve(), args.new.resolve(),
                               args.nruns, args.mpirun)
            (passed if ok else failed).append(test.name)
        except Exception as e:
            print(f"\n  [!] Unexpected error in {test.name}: {e}")
            failed.append(test.name)

    # ── summary ──────────────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("SUMMARY")
    print('='*70)
    if passed:
        print(f"  PASSED ({len(passed)}): {', '.join(passed)}")
    if skipped:
        print(f"  SKIPPED ({len(skipped)}): {', '.join(skipped)}")
    if failed:
        print(f"  FAILED ({len(failed)}): {', '.join(failed)}")
    else:
        print("\n  All completed tests passed.")


if __name__ == "__main__":
    main()
