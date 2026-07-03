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

When both builds include the same floating-point code and only differ by
the bit-identical framework commits, require exact equality instead of
the ULP tolerance, and check the new runtime option:

    python3 scripts/compare_builds.py \\
        --old  /path/to/b8d46fb7/hermes-3 \\
        --new  /path/to/HEAD/hermes-3 \\
        --tests-dir tests/integrated \\
        --expect identical \\
        --verify-check-state-values

Recommended tests for the performance branch changes
-----------------------------------------------------
The performance branch contains two classes of change, with different
expected differences. Choose --expect accordingly:

1. Floating-point reordering (braginskii_collisions, neutral_mixed, ADAS
   cellAverage, integrate.hxx). Same equations, different order of
   operations: differences of ≤ 2-3 ULPs (~1e-15 relative) are expected.
   Use the default --expect ulp when the comparison spans these commits
   (e.g. --old built from master).

2. Bit-identical framework changes (commits 04207481 and 3bb5babf):
   permissions match caching, GuardedOptions session IDs, cached
   field-aligned metrics in div_ops, persistent state Options tree,
   merged guard-cell communication in neutral_full_velocity, and the
   hermes:check_state_values option. These must produce EXACTLY the
   same output, bit for bit, at every timestep. Use --expect identical
   when both binaries include the class-1 commits (e.g. --old built
   from b8d46fb7, --new from HEAD). Any nonzero difference is a bug.

Which tests exercise what:

  Framework changes (state tree, permissions cache, GuardedOptions):
    every test, on every RHS evaluation

  div_ops cached aligned metrics (Div_a_Grad_perp_flows via
  neutral_mixed / classical_diffusion):
    tests/integrated/neutral_mixed
    tests/integrated/collfreq-braginskii-afn
    tests/integrated/collfreq-multispecies
    tests/integrated/2D-production   [requires Zenodo download]
    tests/integrated/2D-recycling    [requires Zenodo download]

  Persistent state tree with the fields:phi / isSection probe paths:
    tests/integrated/vorticity
    tests/integrated/drift-wave
    tests/integrated/alfven-wave

  braginskii_collisions + ADAS reactions (class-1 changes):
    tests/integrated/1D-recycling
    tests/integrated/1D-recycling-dthe

  Multi-processor guard-cell communication (nproc=10; the only MPI
  coverage in the suite, so run at least one of these):
    tests/integrated/2D-production   [requires Zenodo download]
    tests/integrated/2D-recycling    [requires Zenodo download]

  NOT covered by any integrated test: neutral_full_velocity (the merged
  communicate call). Its correctness argument is structural only —
  guard-cell exchanges of independent fields commute.

The hermes:check_state_values option (new in 3bb5babf) can be verified
with --verify-check-state-values: this runs the NEW binary once more
with hermes:check_state_values=false and requires bitwise-identical
output, reporting the time saved by skipping the finiteness sweeps.
Only meaningful for builds with CHECK >= 1.

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
* Differences confined to the x-boundary guard columns of diagnostic fields
  (sources, ddt(...)) are reported but not counted as failures, in BOTH
  --expect modes: Hermes never writes those cells, so the output there is
  stale recycled-allocation content in BOTH builds. Verified on
  2D-production/2D-recycling (the shifted-metric tests, where the div_ops
  metric cache changes the allocation sequence): all interior values and
  all evolved-field values, guard cells included, are bit-identical (or
  within tolerance); only diagnostic x-guard columns hold different
  (junk-scale, e.g. 1e+219) stale values. This applies whenever the two
  binaries being compared differ in allocation sequence for these fields
  (e.g. master vs. any build including 04207481), not just b8d46fb7..HEAD.
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
    # Wall-clock timing — always differ between runs
    'wall_time', 'wtime', 'wtime_comms', 'wtime_io',
    'wtime_per_rhs', 'wtime_per_rhs_e', 'wtime_per_rhs_i',
    'wtime_rhs', 'wtime_invert',
    # RHS call counters — integer counts; differ when the adaptive
    # timestepper takes different numbers of steps (chaotic divergence),
    # not a physics regression
    'ncalls', 'ncalls_e', 'ncalls_i',
    # Solver iteration counter
    'iteration',
}


# Tests with a density/pressure feedback controller are numerically
# sensitive: a bit-level rounding change, even a mathematically valid
# reordering of floating-point operations, can eventually land the density
# trajectory on the opposite side of the controller's threshold. From that
# point on, fields diverge chaotically rather than incorrectly, which
# swamps the diff table with noise unrelated to correctness. For these
# tests, compare only up to the last output timestep observed (in a
# reference run) before that threshold-driven divergence sets in.
_FEEDBACK_SENSITIVE_MAX_TIMESTEP = {
    "1D-recycling": 8,
    "1D-recycling-dthe": 8,
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
    # Also check for restart files, but only for tests that actually fetch
    # them as external input (e.g. from Zenodo). Some runtest scripts just
    # delete any pre-existing BOUT.restart.*.nc as part of their cleanup
    # step, which isn't a signal that the file is a required input.
    runtest = test_dir / "runtest"
    if runtest.is_file():
        text = runtest.read_text()
        if "BOUT.restart" in text and ("zenodo" in text.lower() or "urllib" in text):
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
               mpirun: str,
               extra_args: list[str] | None = None) -> float:
    """Run hermes-3 from work_dir and return wall-clock seconds."""
    if nproc > 1:
        cmd = mpirun.split() + [str(nproc), "./hermes-3", "-d", "data"]
    else:
        cmd = ["./hermes-3", "-d", "data"]
    if extra_args:
        cmd += list(extra_args)

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

def _strip_x_boundaries(da, mxg: int):
    """Drop the x-boundary (guard) columns from a DataArray, if present."""
    if mxg > 0 and "x" in da.dims and da.sizes["x"] > 2 * mxg:
        return da.isel(x=slice(mxg, -mxg))
    return da


def _guard_cell_only(da_old, da_new, mxg: int, ulp_tol: int) -> bool:
    """True if a flagged ULP/rel-diff violation disappears once the
    x-boundary guard columns (never written by Hermes diagnostics, so
    they hold stale recycled-allocation contents in both builds) are
    stripped out. Mirrors the bitwise check in _compare_identical, but
    for the tolerance-based (ULP) comparison path: instead of requiring
    exact equality in the interior, it just requires the interior to
    fall back within the normal tolerance.
    """
    stripped_old = _strip_x_boundaries(da_old, mxg)
    if stripped_old.sizes.get("x", -1) == da_old.sizes.get("x", -2):
        return False  # nothing was actually stripped; not a boundary effect

    ai = numpy.ascontiguousarray(stripped_old.values, dtype=numpy.float64).ravel()
    bi = numpy.ascontiguousarray(_strip_x_boundaries(da_new, mxg).values,
                                 dtype=numpy.float64).ravel()

    mask = numpy.isfinite(ai) & numpy.isfinite(bi) & (ai != 0)
    if not mask.any():
        return True

    a_m, b_m = ai[mask], bi[mask]
    rel = numpy.abs((b_m - a_m) / a_m)
    max_ulps = int(_ulp_diff(a_m, b_m).max())
    return rel.max() <= 1e-10 and max_ulps <= ulp_tol


def _compare_identical(ds_old, ds_new, common) -> bool:
    """Bitwise comparison of every common variable at every timestep.

    Used with --expect identical: the builds under comparison claim to
    perform exactly the same floating-point operations, so the entire
    output trajectory must match bit for bit. Comparing all timesteps
    (rather than just the last one) also makes the feedback-controller
    timestep cap unnecessary — truly identical runs cannot diverge.

    Exception: diffs confined to the x-boundary guard columns are
    reported but do not fail the comparison. Hermes-3 never writes the
    x-guard cells of diagnostic fields (sources, ddt(...)), so their
    output contents are whatever the recycled Array allocation last
    held — e.g. changing how often a temporary is allocated inside an
    operator changes which stale values appear there. They are
    indeterminate in *both* builds and never feed back into the
    simulation (verified: all evolved fields, including their guard
    cells, remain bit-identical).
    """
    nt = min(ds_old.sizes.get('t', 1), ds_new.sizes.get('t', 1))
    if ds_old.sizes.get('t', 1) != ds_new.sizes.get('t', 1):
        print(f"  [!] Different number of output timesteps: "
              f"OLD={ds_old.sizes.get('t', 1)} NEW={ds_new.sizes.get('t', 1)}; "
              f"comparing the first {nt}. This alone indicates the runs "
              f"are not identical.")

    # x-guard columns are only present in the dataset if xhermes kept them
    # (its default); only then can diffs be attributed to them.
    mxg = int(ds_old.metadata.get("MXG", 0)) \
        if ds_old.metadata.get("keep_xboundaries", True) else 0

    print(f"\n  {'Variable':<30} {'differing values':>18}  {'max ULPs':>10}")
    print("  " + "-" * 62)

    any_diff = ds_old.sizes.get('t', 1) != ds_new.sizes.get('t', 1)
    n_identical = 0
    boundary_only = []
    for v in common:
        da_old, da_new = ds_old[v], ds_new[v]
        if 't' in da_old.dims and 't' in da_new.dims:
            da_old = da_old.isel(t=slice(0, nt))
            da_new = da_new.isel(t=slice(0, nt))

        a = numpy.ascontiguousarray(da_old.values, dtype=numpy.float64).ravel()
        b = numpy.ascontiguousarray(da_new.values, dtype=numpy.float64).ravel()

        if a.shape != b.shape:
            print(f"  {v:<30} {'shape mismatch':>18}  <-- DIFFERS")
            any_diff = True
            continue

        # Bitwise comparison catches everything the masked relative-error
        # comparison cannot: differing zeros, signed zeros, NaN payloads
        # and infinities.
        ndiff = int((a.view(numpy.int64) != b.view(numpy.int64)).sum())
        if ndiff == 0:
            n_identical += 1
            continue

        # Are all diffs in the x-boundary guard columns (never written by
        # Hermes diagnostics, so indeterminate in both builds)?
        ai = numpy.ascontiguousarray(_strip_x_boundaries(da_old, mxg).values,
                                     dtype=numpy.float64).ravel()
        bi = numpy.ascontiguousarray(_strip_x_boundaries(da_new, mxg).values,
                                     dtype=numpy.float64).ravel()
        if ai.size < a.size and (ai.view(numpy.int64) == bi.view(numpy.int64)).all():
            boundary_only.append(v)
            print(f"  {v:<30} {ndiff:>18}  {'':>10}  <-- x-guard cells only")
            continue

        any_diff = True
        mask = numpy.isfinite(a) & numpy.isfinite(b)
        max_ulps = int(_ulp_diff(a[mask], b[mask]).max()) if mask.any() else -1
        ulps_str = str(max_ulps) if max_ulps >= 0 else "n/a"
        print(f"  {v:<30} {ndiff:>18}  {ulps_str:>10}  <-- DIFFERS")

    print(f"  ({n_identical}/{len(common)} variables bit-identical over "
          f"all {nt} timestep(s))")
    if boundary_only:
        print(f"  [note] {len(boundary_only)} variable(s) differ only in the "
              f"x-boundary guard columns: {', '.join(boundary_only)}.\n"
              f"         Hermes never writes diagnostic x-guard cells, so both "
              f"builds output stale allocation contents there; not counted "
              f"as a failure.")
    return not any_diff


def compare_outputs(data_old: pathlib.Path, data_new: pathlib.Path,
                    ulp_tol: int, timestep: int = -1,
                    expect_identical: bool = False) -> bool:
    """Compare all field variables at the given timestep. Returns True if ok.

    With expect_identical=True, instead require every common variable to
    be bit-for-bit identical at every timestep.
    """
    ds_old = ds_new = None
    try:
        ds_old = xhermes.open(data_old, unnormalise=False)
        ds_new = xhermes.open(data_new, unnormalise=False)

        if expect_identical:
            vars_old = set(ds_old.data_vars)
            vars_new = set(ds_new.data_vars)
            common = sorted((vars_old & vars_new) - _BOUT_TIMING_VARS)
            sets_differ = bool((vars_old ^ vars_new) - _BOUT_TIMING_VARS)
            if sets_differ:
                print(f"  [!] Variable sets differ: "
                      f"only-OLD={sorted(vars_old - vars_new - _BOUT_TIMING_VARS)} "
                      f"only-NEW={sorted(vars_new - vars_old - _BOUT_TIMING_VARS)}")
            return _compare_identical(ds_old, ds_new, common) and not sets_differ

        nt = min(ds_old.sizes['t'], ds_new.sizes['t'])
        ts = timestep if timestep >= 0 else nt + timestep
        ts = max(0, min(ts, nt - 1))

        old_last = ds_old.isel(t=ts)
        new_last = ds_new.isel(t=ts)

        vars_old = set(old_last.data_vars)
        vars_new = set(new_last.data_vars)
        common = sorted((vars_old & vars_new) - _BOUT_TIMING_VARS)

        if vars_old - vars_new - _BOUT_TIMING_VARS:
            print(f"  [!] Variables only in OLD: "
                  f"{sorted(vars_old - vars_new - _BOUT_TIMING_VARS)}")
        if vars_new - vars_old - _BOUT_TIMING_VARS:
            print(f"  [!] Variables only in NEW: "
                  f"{sorted(vars_new - vars_old - _BOUT_TIMING_VARS)}")

        # x-guard columns are only present if xhermes kept them (its default).
        mxg = int(ds_old.metadata.get("MXG", 0)) \
            if ds_old.metadata.get("keep_xboundaries", True) else 0

        print(f"\n  {'Variable':<30} {'max |rel diff|':>16}  {'RMS rel diff':>14}"
              f"  {'max ULPs':>10}")
        print("  " + "-" * 74)

        any_large_diff = False
        boundary_only = []
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

            if flag and _guard_cell_only(old_last[v], new_last[v], mxg, ulp_tol):
                boundary_only.append(v)
                print(f"  {v:<30} {max_rel:>16.3e}  {rms_rel:>14.3e}"
                      f"  {max_ulps:>10}  <-- x-guard cells only")
                continue

            if flag:
                any_large_diff = True

            print(f"  {v:<30} {max_rel:>16.3e}  {rms_rel:>14.3e}"
                  f"  {max_ulps:>10}{flag}")

        if boundary_only:
            print(f"  [note] {len(boundary_only)} variable(s) differ only in the "
                  f"x-boundary guard columns: {', '.join(boundary_only)}.\n"
                  f"         Hermes never writes diagnostic x-guard cells, so both "
                  f"builds output stale allocation contents there; not counted "
                  f"as a failure.")

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
                  ulp_tol: int,
                  expect_identical: bool = False,
                  verify_check_state_values: bool = False):
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
        if expect_identical:
            # Bit-identical builds cannot diverge, so the whole trajectory
            # is compared and the feedback-controller timestep cap is not
            # needed (a chaotic divergence would itself prove the builds
            # are not identical).
            print(f"\n  Field differences (all timesteps, bitwise):")
            ok = compare_outputs(work_old / "data", work_new / "data", ulp_tol,
                                 expect_identical=True)
            if ok:
                print("\n  => Output is bit-identical.")
            else:
                print("\n  => WARNING: output is NOT bit-identical (see above).")
        else:
            max_ts = _FEEDBACK_SENSITIVE_MAX_TIMESTEP.get(name)
            if max_ts is not None:
                print(f"\n  [!] {name} has a feedback controller that chaotically "
                      f"amplifies any bit-level rounding change late in the run; "
                      f"comparing timestep {max_ts} instead of the last one to "
                      f"avoid that expected (non-correctness) divergence.")
                timestep = max_ts
                label = f"timestep {max_ts}"
            else:
                timestep = -1
                label = "last timestep"
            print(f"\n  Field differences ({label}, ULP tolerance: {ulp_tol}):")
            ok = compare_outputs(work_old / "data", work_new / "data", ulp_tol,
                                 timestep=timestep)
            if ok:
                print("\n  => All variables within tolerance.")
            else:
                print("\n  => WARNING: differences exceed tolerance (see above).")

        # ── hermes:check_state_values verification ───────────────────────────
        if verify_check_state_values:
            print(f"\n  Verifying hermes:check_state_values=false "
                  f"(NEW binary, 1 run)...")
            work_ncsv = tmp / "no_check_state_values"
            setup_work_dir(test_dir, work_ncsv, exec_new)
            try:
                t_ncsv = run_hermes(work_ncsv, nproc, mpirun,
                                    extra_args=["hermes:check_state_values=false"])
            except RuntimeError as e:
                print(f"  [!] Run with check_state_values=false failed: {e}")
                return False
            # Values must be unchanged; only the failure behaviour on
            # non-finite data differs. Compare bitwise against the NEW
            # binary's own (final timed run) output.
            ok_ncsv = compare_outputs(work_new / "data", work_ncsv / "data",
                                      ulp_tol, expect_identical=True)
            saved = mean_new - t_ncsv
            print(f"    check_state_values=false: {t_ncsv:.2f}s vs NEW mean "
                  f"{mean_new:.2f}s ({saved:+.2f}s, "
                  f"{saved / mean_new * 100:+.1f}% of runtime)")
            if ok_ncsv:
                print("    => Identical output with finiteness checks disabled.")
            else:
                print("    => WARNING: disabling check_state_values changed "
                      "the output — this is a bug.")
            ok = ok and ok_ncsv

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
                             "optimisations should stay within 2-3 ULPs. "
                             "Ignored with --expect identical.")
    parser.add_argument("--expect", choices=("ulp", "identical"), default="ulp",
                        help="Expected difference between the builds. "
                             "'ulp' (default): allow float-reordering noise "
                             "up to --ulp-tol at the compared timestep. "
                             "'identical': require bit-for-bit identical "
                             "output at every timestep — use when the "
                             "comparison spans only bit-identical commits "
                             "(e.g. b8d46fb7..HEAD).")
    parser.add_argument("--verify-check-state-values", action="store_true",
                        help="Also run the NEW binary once per test with "
                             "hermes:check_state_values=false, require "
                             "bitwise-identical output, and report the time "
                             "saved. Requires the NEW build to include the "
                             "option (3bb5babf) and CHECK >= 1 for a "
                             "meaningful timing difference.")
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
                expect_identical=(args.expect == "identical"),
                verify_check_state_values=args.verify_check_state_values,
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
