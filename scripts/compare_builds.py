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

Recommended tests for the performance branch changes
-----------------------------------------------------
The following tests directly exercise the modified components:

  braginskii_collisions + braginskii_ion_viscosity + neutral_mixed:
    tests/integrated/collfreq-braginskii-afn
    tests/integrated/collfreq-multispecies

  neutral_mixed (standalone):
    tests/integrated/neutral_mixed

  braginskii_collisions + ADAS reactions (cellAverage/cellAverageInto):
    tests/integrated/1D-recycling
    tests/integrated/1D-recycling-dthe

  neutral_mixed with flux_limit (sqrt reuse path):
    tests/integrated/2D-production   [requires Zenodo download]
    tests/integrated/2D-recycling    [requires Zenodo download]

The optimisations are purely computational (same equations, different order of
floating-point operations), so all field differences should be ≤ 2-3 ULPs
(~1e-15 relative).  Any larger difference indicates a regression.

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
* Use --warmup 1 (the default) to discard the first run, which is slower due
  to cold OS page cache and process start-up overhead.
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

try:
    import xhermes
except ImportError:
    sys.exit("xhermes is required: pip install xhermes")


# Variables BOUT++ writes into every dump file that measure wall-clock timing.
# They will always differ between two separate runs and must not be used for
# correctness comparisons.
_BOUT_TIMING_VARS = {
    'wall_time', 'wtime', 'wtime_comms', 'wtime_io',
    'wtime_per_rhs', 'wtime_per_rhs_e', 'wtime_per_rhs_i',
    'wtime_rhs', 'wtime_invert',
}


# ─── ULP comparison ──────────────────────────────────────────────────────────

def _ulp_diff(a: numpy.ndarray, b: numpy.ndarray) -> numpy.ndarray:
    """Return the number of ULPs between each pair of float64 values.

    Both arrays must contain only finite values.  The algorithm reinterprets
    the IEEE 754 bit pattern as a signed integer, flips the ordering for
    negative values (where the sign-magnitude representation reverses the
    integer ordering), then returns the absolute integer difference.
    """
    INT64_MIN = numpy.int64(-0x8000000000000000)
    ai = a.view(numpy.int64).copy()
    bi = b.view(numpy.int64).copy()
    # Negative floats are ordered backwards in two's-complement; flip them.
    ai[ai < 0] = INT64_MIN - ai[ai < 0]
    bi[bi < 0] = INT64_MIN - bi[bi < 0]
    return numpy.abs(ai - bi)


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


def _clear_dmp_files(data_dir: pathlib.Path) -> None:
    for f in data_dir.glob("BOUT.dmp*.nc"):
        f.unlink()


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

def compare_outputs(data_old: pathlib.Path, data_new: pathlib.Path,
                    ulp_tol: int) -> bool:
    """Compare all field variables in the last timestep. Returns True if ok."""
    ds_old = ds_new = None
    try:
        ds_old = xhermes.open(data_old, unnormalise=False)
        ds_new = xhermes.open(data_new, unnormalise=False)

        old_last = ds_old.isel(t=-1)
        new_last = ds_new.isel(t=-1)

        vars_old = set(old_last.data_vars)
        vars_new = set(new_last.data_vars)
        common = sorted((vars_old & vars_new) - _BOUT_TIMING_VARS)

        if vars_old - vars_new - _BOUT_TIMING_VARS:
            print(f"  [!] Variables only in OLD: "
                  f"{sorted(vars_old - vars_new - _BOUT_TIMING_VARS)}")
        if vars_new - vars_old - _BOUT_TIMING_VARS:
            print(f"  [!] Variables only in NEW: "
                  f"{sorted(vars_new - vars_old - _BOUT_TIMING_VARS)}")

        print(f"\n  {'Variable':<30} {'max |rel diff|':>16}  {'RMS rel diff':>14}"
              f"  {'max ULPs':>10}")
        print("  " + "-" * 74)

        any_large_diff = False
        for v in common:
            a = old_last[v].values.ravel().astype(numpy.float64)
            b = new_last[v].values.ravel().astype(numpy.float64)

            mask = numpy.isfinite(a) & numpy.isfinite(b) & (a != 0)
            if not mask.any():
                continue

            a_m, b_m = a[mask], b[mask]
            rel = numpy.abs((b_m - a_m) / a_m)
            max_rel = rel.max()
            rms_rel = numpy.sqrt((rel**2).mean())
            max_ulps = int(_ulp_diff(a_m, b_m).max())

            rel_flag = "  <-- LARGE REL" if max_rel > 1e-10 else ""
            ulp_flag = f"  <-- HIGH ULPs (>{ulp_tol})" if max_ulps > ulp_tol else ""
            flag = rel_flag or ulp_flag

            if flag:
                any_large_diff = True

            print(f"  {v:<30} {max_rel:>16.3e}  {rms_rel:>14.3e}"
                  f"  {max_ulps:>10}{flag}")

        return not any_large_diff

    except Exception as e:
        print(f"  [!] Could not compare output: {e}")
        return False
    finally:
        for ds in (ds_old, ds_new):
            if ds is not None:
                try:
                    ds.close()
                except Exception:
                    pass


# ─── per-test driver ─────────────────────────────────────────────────────────

_SKIPPED = object()  # sentinel returned when a test is skipped


def run_test_case(test_dir: pathlib.Path,
                  exec_old: pathlib.Path,
                  exec_new: pathlib.Path,
                  nruns: int,
                  warmup: int,
                  mpirun: str,
                  ulp_tol: int):
    """Return True (pass), False (fail), or _SKIPPED."""
    name = test_dir.name
    print(f"\n{'='*70}")
    print(f"Test: {name}  (nproc={detect_nproc(test_dir)})")
    print('='*70)

    if not (test_dir / "data" / "BOUT.inp").is_file():
        print(f"  [!] No data/BOUT.inp — skipping")
        return _SKIPPED

    missing = missing_input_files(test_dir)
    if missing:
        print(f"  [!] Skipping: missing input file(s): {missing}")
        print(f"       Run the original runtest once to download them:")
        print(f"       python3 {test_dir}/runtest")
        return _SKIPPED

    nproc = detect_nproc(test_dir)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        work_old = tmp / "old"
        work_new = tmp / "new"

        setup_work_dir(test_dir, work_old, exec_old)
        setup_work_dir(test_dir, work_new, exec_new)

        def timed_runs(work_dir: pathlib.Path, label: str) -> list[float]:
            times = []
            data_dir = work_dir / "data"

            if warmup > 0:
                print(f"  {label}: warming up ({warmup} run(s))...")
                for _ in range(warmup):
                    try:
                        run_hermes(work_dir, nproc, mpirun)
                    except RuntimeError as e:
                        raise RuntimeError(f"{label} warmup failed: {e}") from e
                    _clear_dmp_files(data_dir)

            for i in range(nruns):
                try:
                    t = run_hermes(work_dir, nproc, mpirun)
                except RuntimeError as e:
                    raise RuntimeError(f"{label} run {i+1} failed: {e}") from e
                times.append(t)
                if i < nruns - 1:
                    _clear_dmp_files(data_dir)
            return times

        # ── old binary ───────────────────────────────────────────────────────
        try:
            times_old = timed_runs(work_old, "OLD")
        except RuntimeError as e:
            print(f"  [!] {e}")
            return False

        # ── new binary ───────────────────────────────────────────────────────
        try:
            times_new = timed_runs(work_new, "NEW")
        except RuntimeError as e:
            print(f"  [!] {e}")
            return False

        # ── timing ───────────────────────────────────────────────────────────
        mean_old = numpy.mean(times_old)
        mean_new = numpy.mean(times_new)
        min_old = numpy.min(times_old)
        min_new = numpy.min(times_new)
        speedup_mean = mean_old / mean_new
        speedup_min = min_old / min_new

        def _sign(x):
            return "faster" if x > 1 else "slower"

        print(f"\n  Timing ({nruns} timed run(s) each, {warmup} warmup):")
        print(f"    OLD: mean={mean_old:.2f}s  min={min_old:.2f}s"
              f"  ({', '.join(f'{t:.2f}' for t in times_old)})")
        print(f"    NEW: mean={mean_new:.2f}s  min={min_new:.2f}s"
              f"  ({', '.join(f'{t:.2f}' for t in times_new)})")
        print(f"    Speedup (mean): {speedup_mean:.3f}x"
              f" ({abs(speedup_mean-1)*100:.1f}% {_sign(speedup_mean)})")
        print(f"    Speedup (min):  {speedup_min:.3f}x"
              f" ({abs(speedup_min-1)*100:.1f}% {_sign(speedup_min)})")

        # ── correctness ──────────────────────────────────────────────────────
        print(f"\n  Field differences (last timestep, ULP tolerance: {ulp_tol}):")
        ok = compare_outputs(work_old / "data", work_new / "data", ulp_tol)
        if ok:
            print("\n  => All variables within tolerance.")
        else:
            print("\n  => WARNING: differences exceed tolerance (see above).")

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
    parser.add_argument("--warmup", type=int, default=1,
                        help="Throwaway runs before timing to warm the page "
                             "cache and stabilise clocks (default: 1)")
    parser.add_argument("--ulp-tol", type=int, default=10,
                        help="Flag variables whose max ULP difference exceeds "
                             "this threshold (default: 10). Pure computational "
                             "optimisations should stay within 2-3 ULPs.")
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
            result = run_test_case(
                test,
                args.old.resolve(),
                args.new.resolve(),
                nruns=args.nruns,
                warmup=args.warmup,
                mpirun=args.mpirun,
                ulp_tol=args.ulp_tol,
            )
            if result is _SKIPPED:
                skipped.append(test.name)
            elif result:
                passed.append(test.name)
            else:
                failed.append(test.name)
        except Exception as e:
            print(f"\n  [!] Unexpected error in {test.name}: {e}")
            failed.append(test.name)

    # ── summary ──────────────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("SUMMARY")
    print('='*70)
    if passed:
        print(f"  PASSED  ({len(passed)}): {', '.join(passed)}")
    if skipped:
        print(f"  SKIPPED ({len(skipped)}): {', '.join(skipped)}")
    if failed:
        print(f"  FAILED  ({len(failed)}): {', '.join(failed)}")
    if not failed:
        print("\n  All completed tests passed.")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
