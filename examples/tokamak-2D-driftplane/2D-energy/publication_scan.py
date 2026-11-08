#!/usr/bin/env python3
"""Matched interchange resolution study. See publication-study.md for usage."""

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import signal
import string
import subprocess
import sys
import tempfile
import time
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(os.environ.get("TMPDIR", "/tmp")) / "hermes-publication-matplotlib"),
)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE_ROOT = ROOT
DEFAULT_EXE = (
    ROOT / "hermes-3" if (ROOT / "hermes-3").is_file() else ROOT / "build/hermes-3"
)
sys.path.insert(0, str(ROOT / "tools"))
from energy_analysis import (
    AuditError,
    PrescribedSources,
    energy_densities,
    energy_inventory,
    open_run,
    regular_timestep,
)
from energy_analysis.convergence import error_norms, global_differences, restrict
from energy_analysis.fluxes import fill_wall_guards, radial_boundary_mode
from energy_analysis.operators import perpendicular_operator
from energy_analysis.pressure import legacy_operator, legacy_viscous_work
from energy_analysis.report import revision

SCHEMA = 1
DIFFERENCES = [
    "Pressure uses a potential-flow polarisation approximation and electron compression; total energy uses ion-current partition and species transport.",
    "Pressure polarisation potential includes electron pressure with its signed mass/charge weight; total energy includes ions only.",
    "Pressure native inversion uses FD X / spectral Z; total energy uses face/corner FV X-Z.",
    "Pressure retains cell-centred corner-gradient viscous heating; total energy uses conservative face work.",
    "The pressure perturbation vanishes at both X faces. Pressure uses Neumann ion-pressure guards; total energy imposes mean-Neumann, zero-AC-midpoint pressure and psi guards. Physical initial fields are matched and boundary-adjacent cells remain in the error norms.",
    "Explicit SI metric pairs enforce the 5 cm box; the legacy recalculate_metric path produces inconsistent Z metric pairs for Bnorm=0.35 T.",
    "Common phi solves the continuum ion-only zero-vorticity problem; both discrete vorticities are initialized from this same sampled physical state.",
]


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def runtime_artifacts(exe):
    """Hash the executable and any separately built Hermes/BOUT libraries."""
    artifacts = {exe.resolve()}
    for directory in (
        exe.parent,
        exe.parent / "lib",
        exe.parent / "external/BOUT-dev/lib",
    ):
        for prefix in ("libhermes-3-lib", "libbout++"):
            for path in directory.glob(prefix + "*"):
                if path.is_file() and (
                    path.suffix == ".dylib"
                    or ".so" in path.name
                    or path.suffix == ".dll"
                ):
                    artifacts.add(path.resolve())
    return {str(path): digest(path) for path in sorted(artifacts)}


def save_json(path, value):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def write_netcdf(dataset, path):
    # xBOUT attaches Python metadata/options objects to coordinates as well as
    # fields. Keep only plain netCDF attributes in retained analysis arrays.
    dataset = dataset.copy(deep=False)
    simple = (str, int, float, np.number)
    dataset.attrs = {k: v for k, v in dataset.attrs.items() if isinstance(v, simple)}
    for name in dataset.variables:
        dataset[name].attrs = {
            k: v for k, v in dataset[name].attrs.items() if isinstance(v, simple)
        }
    dataset.to_netcdf(path)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    p.add_argument("--levels", type=int, default=4)
    p.add_argument("--base-nx", type=int, default=32)
    p.add_argument("--base-nz", type=int, default=32)
    p.add_argument("--refinement", type=int, default=2)
    p.add_argument(
        "--phi-boundary", choices=("fixed", "relax", "both"), default="fixed"
    )
    p.add_argument(
        "--phi-boundary-timescale",
        type=float,
        default=1e-4,
        help="Legacy vorticity relaxation timescale [s]",
    )
    p.add_argument(
        "--phi-boundary-velocity",
        type=float,
        default=1000,
        help="Total-ion-energy relaxation velocity [m/s]",
    )
    p.add_argument("--viscosity", type=float, default=1)
    p.add_argument("--end-time", type=float, default=20000)
    p.add_argument("--outputs", type=int, default=200)
    p.add_argument("--snapshot-times", type=float, nargs=4)
    p.add_argument("--convergence-time", type=float)
    p.add_argument(
        "--field",
        choices=("ion_pressure", "vorticity", "density"),
        default="ion_pressure",
    )
    p.add_argument("--nproc", type=int, default=1)
    p.add_argument("--timeout", type=float, default=7200)
    p.add_argument("--output-root", type=Path, default=HERE / "publication-results")
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--analyse-only", type=Path, metavar="STUDY")
    mode.add_argument("--resume", type=Path, metavar="STUDY")
    mode.add_argument(
        "--extend",
        type=Path,
        metavar="STUDY",
        help="Inherit an existing study's settings and append grids using --levels",
    )
    mode.add_argument(
        "--import-runs",
        type=Path,
        metavar="STUDY",
        help="Validate returned outputs and register them without running Hermes",
    )
    p.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare matched time-zero restarts and a standalone remote runner, then stop",
    )
    p.add_argument(
        "--cases",
        nargs="+",
        metavar="CASE",
        help="Select case directory names for preparation, execution or import",
    )
    p.add_argument(
        "--import-from",
        type=Path,
        metavar="DIRECTORY",
        help="Parent directory of returned case directories; defaults to STUDY",
    )
    p.add_argument(
        "--allow-runner-change",
        action="store_true",
        help="With --resume, record a changed runner version; input and executable checks still apply",
    )
    p.add_argument("--diagnose", action="store_true")
    p.add_argument("--atol", type=float, default=1e-10)
    p.add_argument("--rtol", type=float, default=1e-5)
    p.add_argument("--phi-rtol", type=float, default=1e-5)
    p.add_argument(
        "--phi-atol", type=float, default=1e-3, help="Potential iteration scale [eV]"
    )
    p.add_argument("--vort-atol", type=float, default=1e-10)
    p.add_argument("--max-timestep", type=float, default=50)
    p.add_argument(
        "--accuracy-pilot",
        action="store_true",
        help="Repeat two finest grids with tightened tolerances/timesteps and output cadence",
    )
    p.add_argument("--latex", action="store_true")
    p.add_argument("--style", type=Path, help="JSON Matplotlib rcParams overrides")
    p.add_argument(
        "--energy-format",
        choices=("png", "pdf-panels"),
        default="png",
        help="Energy histories: one six-panel PNG (default) or six standalone PDFs",
    )
    p.add_argument(
        "--smoke",
        action="store_true",
        help="Two 8/16 grids, both boundaries, inviscid, end time 12, 4 intervals",
    )
    return p


def configuration(args):
    if args.smoke:
        args.levels, args.base_nx, args.base_nz = 2, 8, 8
        args.phi_boundary, args.viscosity, args.end_time, args.outputs = (
            "both",
            0,
            12,
            4,
        )
    c = {
        k: v
        for k, v in vars(args).items()
        if k
        not in (
            "exe",
            "output_root",
            "resume",
            "extend",
            "import_runs",
            "prepare_only",
            "cases",
            "import_from",
            "allow_runner_change",
            "analyse_only",
            "style",
            "latex",
            "energy_format",
            "smoke",
            "timeout",
        )
    }
    if c["levels"] < 2 or min(c["base_nx"], c["base_nz"]) < 4 or c["refinement"] < 2:
        raise AuditError("Use >=2 grids, >=4 physical cells, refinement>=2")
    for key in (
        "end_time",
        "phi_boundary_timescale",
        "phi_boundary_velocity",
        "atol",
        "rtol",
        "phi_atol",
        "phi_rtol",
        "vort_atol",
        "max_timestep",
    ):
        if not np.isfinite(c[key]) or c[key] <= 0:
            raise AuditError(f"{key} must be finite and positive")
    if (
        not np.isfinite(c["viscosity"])
        or c["viscosity"] < 0
        or c["nproc"] < 1
        or c["outputs"] < 4
    ):
        raise AuditError("Invalid viscosity, process count or output intervals")
    if c["base_nx"] % c["nproc"] or c["base_nx"] // c["nproc"] < 4:
        raise AuditError(
            "X decomposition needs >=4 physical cells per rank on every grid"
        )
    c["snapshot_times"] = c["snapshot_times"] or [
        0,
        c["end_time"] / 4,
        c["end_time"] / 2,
        c["end_time"],
    ]
    c["convergence_time"] = (
        c["convergence_time"]
        if c["convergence_time"] is not None
        else c["end_time"] / 2
    )
    # Refine the regular cadence to include all requested times exactly.
    from fractions import Fraction

    denominator = c["outputs"]
    for value in [*c["snapshot_times"], c["convergence_time"]]:
        if not np.isfinite(value) or not 0 <= value <= c["end_time"]:
            raise AuditError("Requested times must lie in [0,end_time]")
        f = Fraction(value / c["end_time"]).limit_denominator(10000)
        if not np.isclose(float(f), value / c["end_time"], atol=1e-12, rtol=0):
            raise AuditError(
                "Requested times cannot be represented by a bounded regular cadence"
            )
        denominator = math.lcm(denominator, f.denominator)
    if denominator > 10000:
        raise AuditError("Requested snapshot times require >10000 output intervals")
    c["outputs"] = denominator
    return c


def input_text(config, row):
    pressure = row["formulation"] == "pressure"
    relaxing = row["boundary"] == "relax"
    c = config
    # The pressure component uses legacy flags; total ion energy owns its
    # solver boundary flags through radial_boundary. Generate only the options
    # each component accepts, including during time-zero startup.
    if pressure:
        boundary_options = f"phi_boundary_relax = {str(relaxing).lower()}"
        if relaxing:
            boundary_options += (
                f"\nphi_boundary_timescale = {c['phi_boundary_timescale']}"
                "\nphi_core_averagey = false"
            )
    else:
        boundary_options = f"radial_boundary = {'no_slip' if relaxing else 'fixed'}"
        if relaxing:
            boundary_options += (
                f"\nphi_boundary_velocity = {c['phi_boundary_velocity']}"
            )
    values = {
        "mesh_nx": row["nx"] + 4,
        "nz": row["nz"],
        "nproc": c["nproc"],
        "outputs": c["outputs"],
        "output_step": c["end_time"] / c["outputs"],
        "atol": c["atol"],
        "rtol": c["rtol"],
        "max_timestep": c["max_timestep"],
        "viscosity": c["viscosity"],
        "diagnose": str(c["diagnose"]).lower(),
        "boundary_options": boundary_options,
        "laplacian_boundary_options": (
            "inner_boundary_flags = 16\nouter_boundary_flags = 16" if pressure else ""
        ),
        "model": "vorticity" if pressure else "evolve_vorticity_ion_energy",
        "components": "vorticity, polarisation_drift"
        if pressure
        else "evolve_vorticity_ion_energy",
        "ion_components": "evolve_density, evolve_pressure"
        if pressure
        else "evolve_density",
        "ion_pressure_options": "thermal_conduction = false\nlow_p_diffuse_perp = false\nlow_T_diffuse_perp = false\ndamp_p_nt = false\np_div_v = false"
        if pressure
        else "",
        "model_options": (
            "diamagnetic = true\ndiamagnetic_polarisation = true\nbndry_flux = false\nexb_advection = true\nexb_advection_simplified = false\nsplit_n0 = false\nphi_dissipation = false\nvort_dissipation = false\nsheath_boundary = false\nhyper_z = -1"
            if pressure
            else f"density_floor = 1e-7\nphi_atol = {c['phi_atol']}\nphi_rtol = {c['phi_rtol']}\nvort_atol = {c['vort_atol']}\nphi_max_iterations = 400"
        ),
        "polarisation_options": (
            "[polarisation_drift]\nboussinesq = true\naverage_atomic_mass = 1\ndiamagnetic_polarisation = true\nadvection = true\ndiagnose = "
            + str(c["diagnose"]).lower()
        )
        if pressure
        else "",
    }
    text = string.Template(
        (HERE / "publication-input.inp.template").read_text()
    ).substitute(values)
    # Only the formulation's evolving fields have initial-profile sections.
    import re

    section = "Eion" if pressure else "Ph+"
    text = re.sub(
        r"\[" + re.escape(section) + r"\]\n.*?(?=\n\[)", "", text, flags=re.DOTALL
    )
    return text


def check_output_times(values, config, source, *, seed=False):
    expected = (
        np.asarray([0.0])
        if seed
        else np.linspace(0, config["end_time"], config["outputs"] + 1)
    )
    if values.shape != expected.shape or not np.allclose(
        values, expected, rtol=0, atol=max(1, config["end_time"]) * 1e-10
    ):
        raise AuditError(f"Incomplete/mismatched output times: {source}")


def verify_completed_outputs(directory, config, *, seed=False):
    """Check every rank's time axis and final physical state after a bad exit."""
    expected = {f"BOUT.dmp.{rank}.nc" for rank in range(config["nproc"])}
    files = sorted(directory.glob("BOUT.dmp.*.nc"))
    if {file.name for file in files} != expected:
        raise AuditError("Expected one dump file per MPI rank")
    for file in files:
        with xr.open_dataset(file, decode_times=False) as ds:
            if "t" not in ds:
                raise AuditError(f"Missing output time axis: {file}")
            check_output_times(np.asarray(ds.t), config, file, seed=seed)
            fields = ["Ne", "Nh+", "Pe", "phi", "Vort"]
            fields += ["Eion", "Tion"] if "Eion" in ds else ["Ph+"]
            for name in fields:
                if name not in ds or not {"t", "x", "y", "z"}.issubset(ds[name].dims):
                    raise AuditError(
                        f"Missing time-dependent state field {name}: {file}"
                    )
                final = ds[name].isel(t=-1, x=slice(2, -2)).load()
                if not final.size or not np.isfinite(final).all():
                    raise AuditError(f"Incomplete/nonfinite final {name}: {file}")
    if seed:
        restarts = sorted(directory.glob("BOUT.restart.*.nc"))
        if {file.name for file in restarts} != {
            f"BOUT.restart.{rank}.nc" for rank in range(config["nproc"])
        }:
            raise AuditError("Startup did not write one restart per MPI rank")
        for file in restarts:
            with xr.open_dataset(file, decode_cf=False) as ds:
                ds.load()
    return {
        "ranks": config["nproc"],
        "records_per_rank": 1 if seed else config["outputs"] + 1,
        "final_time_normalised": 0.0 if seed else config["end_time"],
    }


def accept_completed_outputs(directory, config, status, *, seed=False):
    path = directory / ("startup-status.json" if seed else "status.json")
    try:
        check = verify_completed_outputs(directory, config, seed=seed)
    except (
        AuditError,
        OSError,
        ValueError,
        KeyError,
        RuntimeError,
        TypeError,
    ) as error:
        status["outputs_complete"] = False
        status["output_check_error"] = str(error)
        save_json(path, status)
        return False
    status["outputs_complete"] = True
    status["output_check"] = check
    status.pop("output_check_error", None)
    save_json(path, status)
    print(
        f"Warning: exit code {status['exit_code']} from {directory}; "
        "all expected outputs were verified. Accepting the completed run.",
        file=sys.stderr,
        flush=True,
    )
    return True


def completed_run(directory, config):
    path = directory / "status.json"
    if not path.exists():
        return False
    status = json.loads(path.read_text())
    if status.get("origin") == "external":
        verify_completed_outputs(directory, config)
        return status.get("outputs_complete", False)
    if status["exit_code"] == 0:
        return True
    # Recheck files rather than trusting a saved completion flag on resume.
    return accept_completed_outputs(directory, config, status)


def execute(exe, directory, config, timeout, *, seed=False):
    command = [str(exe), "-d", str(directory)]
    if seed:
        command += ["restart=false", "solver:nout=0"]
    if config["nproc"] > 1:
        mpi = shutil.which("mpiexec")
        if mpi is None:
            raise AuditError("mpiexec is required for nproc>1")
        command = [mpi, "-n", str(config["nproc"]), *command]
    log = directory / ("startup.log" if seed else "run.log")
    start = time.monotonic()
    interrupted = False
    termination_reason = None
    with log.open("w") as stream:
        process = subprocess.Popen(
            command,
            stdout=stream,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            env={
                **os.environ,
                "FI_SOCKETS_IFACE": os.environ.get("FI_SOCKETS_IFACE", "en0"),
            },
        )
        try:
            code = process.wait(timeout=timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
            interrupted = isinstance(error, KeyboardInterrupt)
            termination_reason = "interrupted" if interrupted else "timeout"
            code = -2 if interrupted else -1
            # Stop the entire local launch group, including MPI children.
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
    status = {
        "command": command,
        "exit_code": code,
        "elapsed_s": time.monotonic() - start,
        "log": log.name,
    }
    if termination_reason:
        status["termination_reason"] = termination_reason
    save_json(directory / ("startup-status.json" if seed else "status.json"), status)
    if interrupted:
        raise KeyboardInterrupt
    if code and not accept_completed_outputs(directory, config, status, seed=seed):
        raise AuditError(
            f"Run failed ({code}): {log}\n"
            + "\n".join(log.read_text().splitlines()[-30:])
        )
    return status


def initial_state(run):
    """Sample the analytic ion-only continuum solution with separate guards."""
    ds = run.data
    nx, nz = ds.sizes["x"] - 4, ds.sizes["z"]
    x = (np.arange(nx + 4) - 2 + 0.5) / nx
    z = np.arange(nz) * 2 * np.pi / nz
    shape = (nx + 4, ds.sizes["y"], nz)
    pressure = 1 + 0.1 * np.sin(2 * np.pi * x[:, None]) * np.sin(
        2 * np.pi * x[:, None] - z[None, :]
    )
    # Both pressure wall values are constant, so their harmonic extension is 1.
    potential = 1 - pressure
    p = xr.DataArray(
        np.broadcast_to(pressure[:, None, :], shape).copy(),
        dims=("x", "y", "z"),
        coords={d: ds[d] for d in ("x", "y", "z")},
    )
    phi = xr.DataArray(
        np.broadcast_to(potential[:, None, :], shape).copy(),
        dims=p.dims,
        coords=p.coords,
    )
    for i in range(2):
        p[{"x": i}] = p.isel(x=3 - i)
        p[{"x": -i - 1}] = p.isel(x=-4 + i)
    phi[{"x": 1}] = -phi.isel(x=2)
    phi[{"x": -2}] = -phi.isel(x=-3)
    phi[{"x": 0}] = phi.isel(x=1)
    phi[{"x": -1}] = phi.isel(x=-2)
    ion_pressure = (
        p if run.formulation == "pressure" else fill_wall_guards(p, run, "pressure")
    )
    psi = p * (1 - 1 / 1836) if run.formulation == "pressure" else ion_pressure
    c = run.scalar_option("average_atomic_mass") / ds.Bxy**2
    op = (
        legacy_operator(c, phi + psi, ds)
        if run.formulation == "pressure"
        else perpendicular_operator(c, phi + psi, ds)
    )
    vort = op.divergence.copy()
    for i in range(2):
        vort[{"x": i}] = vort.isel(x=3 - i)
        vort[{"x": -i - 1}] = vort.isel(x=-4 + i)
    if run.formulation == "total_energy":
        vort = fill_wall_guards(vort, run, "neumann")
    energy = 1.5 * ion_pressure + op.kperp
    for i in range(2):
        energy[{"x": i}] = energy.isel(x=3 - i)
        energy[{"x": -i - 1}] = energy.isel(x=-4 + i)
    return {
        "Pe": p,
        "Ph+": ion_pressure,
        "Nh+": xr.ones_like(p),
        "phi": phi,
        "Vort": vort,
        "Eion": energy,
    }


def build_startup(exe, directory, config, timeout):
    """Obtain native restart metadata, then install matched time-zero arrays."""
    seed = directory / "startup"
    seed.mkdir()
    shutil.copyfile(directory / "BOUT.inp", seed / "BOUT.inp")
    execute(exe, seed, config, timeout, seed=True)
    run = open_run(seed)
    initial = initial_state(replace(run, data=run.data.isel(t=0, drop=True)))
    # Every nonperiodic total-energy run now stores phi midpoints in Vort guards.
    relaxing = run.formulation == "total_energy" or run.data.options[
        run.model
    ].get_bool("phi_boundary_relax")
    nx = run.data.sizes["x"] - 4
    files = sorted(
        seed.glob("BOUT.restart.*.nc"), key=lambda f: int(f.name.split(".")[-2])
    )
    if len(files) != config["nproc"]:
        raise AuditError("Startup did not write one restart per MPI rank")
    for rank, file in enumerate(files):
        with xr.open_dataset(file, decode_cf=False) as opened:
            restart = opened.load()
        for name, field in initial.items():
            if name not in restart:
                continue
            values = np.asarray(field)[
                rank * (nx // config["nproc"]) : rank * (nx // config["nproc"])
                + nx // config["nproc"]
                + 4
            ].copy()
            if name == "Vort" and relaxing:
                if rank == 0:
                    values[:2] = 0
                if rank == config["nproc"] - 1:
                    values[-2:] = 0
            restart[name].data = values
        restart["tt"].data = np.asarray(0.0)
        restart["hist_hi"].data = np.asarray(0, dtype=restart.hist_hi.dtype)
        # BOUT rewrites restart metadata in place during postInit. xarray's
        # automatic NaN fill attributes are not part of the native restart
        # schema and cannot be redefined after data have been written.
        for name in restart.variables:
            restart[name].attrs.pop("_FillValue", None)
        restart.to_netcdf(
            directory / file.name,
            encoding={name: {"_FillValue": None} for name in restart.variables},
        )
    # Reference arrays retained independently from the solver's restart encoding.
    write_netcdf(xr.Dataset(initial), directory / "initial-state.nc")
    run.data.close()


def validate_times(run, config):
    values = np.asarray(run.data.t)
    regular_timestep(values)
    check_output_times(values, config, run.source)
    return [exact_index(values, v) for v in config["snapshot_times"]], exact_index(
        values, config["convergence_time"]
    )


def exact_index(times, value):
    matches = np.flatnonzero(
        np.isclose(times, value, rtol=0, atol=max(1, np.max(abs(times))) * 1e-10)
    )
    if len(matches) != 1:
        raise AuditError(f"Requested time {value} is not an exact unique saved record")
    return int(matches[0])


def field_array(run, index, field):
    ds = run.data
    if field == "perpendicular_kinetic":
        snapshot = replace(run, data=ds.isel(t=[index]).compute())
        value = energy_densities(snapshot).ion_perpendicular.isel(t=0)
        return np.asarray(
            run.interior(value).squeeze("y").transpose("x", "z")
        ) * run.energy_scale, "J m$^{-3}$"
    if field == "ion_pressure":
        value = ds["Ph+"] if run.formulation == "pressure" else ds["Nh+"] * ds.Tion
        scale = run.energy_scale
        unit = "Pa"
    elif field == "density":
        value, scale, unit = ds["Nh+"], ds.metadata["Nnorm"], "m$^{-3}$"
    else:
        from scipy.constants import elementary_charge

        value, scale, unit = (
            ds.Vort,
            elementary_charge * ds.metadata["Nnorm"],
            "C m$^{-3}$",
        )
    return np.asarray(
        run.interior(value.isel(t=index)).squeeze("y").transpose("x", "z")
    ) * scale, unit


def coordinates(run):
    # Derive physical centers from actual normalised metric cell spacings and
    # the BOUT center convention used by the initial profile (Z begins at zero).
    ds = run.data
    nx, nz = ds.sizes["x"] - 4, ds.sizes["z"]
    dx = (
        float(ds.dx.values.flat[0])
        / np.sqrt(float(ds.g11.values.flat[0]))
        * run.length_scale
    )
    dz = (
        float(ds.dz.values.flat[0])
        / np.sqrt(float(ds.g33.values.flat[0]))
        * run.length_scale
    )
    lx, lz = dx * nx, dz * nz
    if not np.allclose([lx, lz], [0.05, 0.05], rtol=1e-10):
        raise AuditError("Physical box changed with resolution")
    return (np.arange(nx) + 0.5) * dx, np.arange(nz) * dz, lx, lz


def validate_initial_state(directory, run, *, first=None):
    """Check the physical time-zero state independently of execution status."""
    with xr.open_dataset(directory / "initial-state.nc") as saved:
        initial = saved.load()
    checks = {}
    if first is None:
        first = replace(run, data=run.data.isel(t=[0]).compute())
    names = ["Pe", "Nh+", "phi", "Vort"]
    if run.formulation == "total_energy":
        names.append("Eion")
    for name in names:
        checks[name] = float(
            abs(run.interior(first.data[name].isel(t=0) - initial[name])).max(
                skipna=False
            )
        )
        if not np.isfinite(checks[name]) or checks[name] > 1e-8:
            raise AuditError(f"Initial {name} mismatch {checks[name]}: {directory}")
    pi = (
        first.data["Ph+"]
        if run.formulation == "pressure"
        else first.data["Nh+"] * first.data.Tion
    )
    checks["ion_pressure"] = float(
        abs(run.interior(pi.isel(t=0) - initial["Ph+"])).max(skipna=False)
    )
    if not np.isfinite(checks["ion_pressure"]) or checks["ion_pressure"] > 1e-8:
        raise AuditError("Initial ion pressure mismatch")
    checks["Ne"] = float(
        abs(run.interior(first.data.Ne.isel(t=0) - initial["Nh+"])).max(skipna=False)
    )
    if not np.isfinite(checks["Ne"]) or checks["Ne"] > 1e-8:
        raise AuditError("Initial electron density mismatch")
    return checks


def analyse_run(directory, config):
    if not completed_run(directory, config):
        raise AuditError(f"Refusing failed run {directory}")
    run = open_run(directory)
    save_json(directory / "resolved-options.json", run.data.options.as_dict())
    snapshots, convergence = validate_times(run, config)
    first = replace(run, data=run.data.isel(t=[0]).compute())
    checks = validate_initial_state(directory, run, first=first)
    inventories, powers, defects, diagnostics = [], [], [], []
    coverage = (
        "Legacy polarisation and full Phi boundary-state power coverage is incomplete"
    )
    boundary_midpoints, boundary_rates = [], []
    sources = PrescribedSources(
        *({s.name: 0 for s in run.species} for _ in range(3)),
        provenance="publication template: all prescribed sources disabled",
    )
    for i in range(run.data.sizes["t"]):
        snapshot = replace(run, data=run.data.isel(t=[i]).compute())
        densities = energy_densities(snapshot)
        inv = energy_inventory(snapshot, densities=densities).compute()
        if not all(np.isfinite(inv[v]).all() for v in inv):
            raise AuditError("Nonfinite energy inventory")
        inventories.append(inv)
        if run.formulation == "total_energy":
            from energy_analysis.fluxes import (
                closed_slab_outward_power,
                compare_flux_diagnostics,
                reconstruct_fluxes,
            )

            boundary = closed_slab_outward_power(
                snapshot, densities=densities
            ).compute()
            power = float(boundary.outward.isel(t=0))
            if not np.isfinite(power):
                raise AuditError("Nonfinite reconstructed boundary power")
            powers.append(power)
            if config["diagnose"]:
                fluxes = reconstruct_fluxes(snapshot, sources=sources).compute()
                diagnostics.append(compare_flux_diagnostics(snapshot, fluxes))
            else:
                diagnostics.append({})
            defects.append(None)
            coverage = "Complete native total-ion plus electron thermal balance; reconstructed radial faces; zero prescribed sources"
        else:
            viscous = legacy_viscous_work(snapshot)
            defects.append(
                float(run.integrate(viscous.viscous_heating_defect).isel(t=0))
                * run.energy_scale
                * run.length_scale**3
                * run.frequency
            )
            if config["diagnose"]:
                from energy_analysis.fluxes import compare_flux_diagnostics

                diagnostics.append(compare_flux_diagnostics(snapshot, viscous))
            else:
                diagnostics.append({})
        phi = snapshot.data.phi.isel(t=0)
        mid = [(phi.isel(x=c) + phi.isel(x=g)) / 2 for c, g in ((2, 1), (-3, -2))]
        if any(float(abs(v - v.mean("z")).max()) > 1e-9 for v in mid):
            raise AuditError("Phi midpoint varies in Z")
        b = [float(v.mean()) for v in mid]
        boundary_midpoints.append(b)
        target = [float(phi.isel(x=c).mean()) for c in (2, -3)]
        if run.formulation == "total_energy":
            if str(run.data.options[run.model]["radial_boundary"]) == "no_slip":
                from energy_analysis.fluxes import boundary_relaxation

                rates = boundary_relaxation(snapshot).phi_boundary_rate
                boundary_rates.append(
                    [
                        float(rates.sel(boundary=side).mean())
                        for side in ("inner", "outer")
                    ]
                )
            else:
                boundary_rates.append([0.0, 0.0])
        else:
            boundary_rates.append(
                [
                    (t - v) / (config["phi_boundary_timescale"] * run.frequency)
                    for t, v in zip(target, b)
                ]
            )
    inventory = xr.concat(inventories, "t")
    write_netcdf(inventory, directory / "inventory.nc")
    native, common = (
        np.asarray(inventory.plasma_native),
        np.asarray(inventory.plasma_reconstructed),
    )
    t = np.asarray(run.data.t) / run.frequency
    cumulative = (
        np.r_[
            0,
            np.cumsum(
                (np.asarray(powers)[1:] + np.asarray(powers)[:-1]) * np.diff(t) / 2
            ),
        ]
        if powers
        else None
    )
    residual = native - native[0] + cumulative if cumulative is not None else None
    drift, common_drift = native - native[0], common - common[0]
    peak = float(np.ptp(np.asarray(inventory.ion_perpendicular)))
    bx = np.asarray(boundary_midpoints)
    if radial_boundary_mode(run) == "fixed" and float(abs(bx).max()) > 1e-9:
        raise AuditError("Fixed phi boundary midpoint changed")
    ode = (
        np.diff(bx, axis=0)
        - (np.asarray(boundary_rates)[1:] + np.asarray(boundary_rates)[:-1])
        * np.diff(np.asarray(run.data.t))[:, None]
        / 2
    )
    partition = (
        float(
            abs(
                run.interior(
                    energy_densities(first).ion_reconstructed - first.data.Eion
                )
            ).max()
        )
        if run.formulation == "total_energy"
        else None
    )
    report = {
        "initial_checks_normalised": checks,
        "coverage": coverage,
        "complete": powers != [],
        "initial_native_J": float(native[0]),
        "final_native_J": float(native[-1]),
        "final_drift_J": float(drift[-1]),
        "max_relative_drift": float(abs(drift).max() / abs(native[0])),
        "common_final_drift_J": float(common_drift[-1]),
        "max_native_common_difference_J": float(abs(native - common).max()),
        "drift_relative_peak_kinetic_exchange": float(abs(drift).max() / peak)
        if peak
        else None,
        "integrated_outward_J": float(cumulative[-1])
        if cumulative is not None
        else None,
        "final_corrected_residual_J": float(residual[-1])
        if residual is not None
        else None,
        "max_corrected_residual_J": float(abs(residual).max())
        if residual is not None
        else None,
        "partition_max_error_normalised": partition,
        "legacy_viscous_defect_W": defects,
        "outward_power_W": powers if powers else None,
        "corrected_residual_J": residual.tolist() if residual is not None else None,
        "boundary_phi_normalised": boundary_midpoints,
        "boundary_ode_interval_error_normalised": float(abs(ode).max()),
        "diagnostics": diagnostics,
        "time_normalised": np.asarray(run.data.t).tolist(),
        "volume_m3": float(run.integrate(xr.ones_like(first.data.phi)).isel(t=0))
        * run.length_scale**3,
    }
    save_json(directory / "analysis.json", report)
    x, z, lx, lz = coordinates(run)
    fields = [field_array(run, i, config["field"])[0] for i in snapshots]
    kinetic_fields = (
        [field_array(run, i, "perpendicular_kinetic")[0] for i in snapshots]
        if run.formulation == "total_energy"
        else []
    )
    comparison, unit = field_array(run, convergence, config["field"])
    volume = np.asarray(
        run.interior(run.data.J * run.data.dx * run.data.dy * run.data.dz).squeeze("y")
    )
    if volume.ndim == 1:
        volume = np.broadcast_to(volume[:, None], comparison.shape)
    q = float(inventory.ion_perpendicular.isel(t=convergence))
    run.data.close()
    return {
        "report": report,
        "inventory": inventory,
        "x": x,
        "z": z,
        "lx": lx,
        "lz": lz,
        "snapshots": fields,
        "kinetic_snapshots": kinetic_fields,
        "comparison": comparison,
        "volume": volume,
        "unit": unit,
        "Q": q,
    }


def write_tables(study, energies, convergence, config, *, latex=False):
    for name, rows in (("energy", energies), ("convergence", convergence)):
        save_json(study / (name + "-table.json"), rows)
        with (study / (name + "-table.csv")).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        fields = (
            [
                "boundary",
                "formulation",
                "grid",
                "final_drift_J",
                "max_relative_drift",
                "final_corrected_residual_J",
            ]
            if name == "energy"
            else [
                "boundary",
                "formulation",
                "grid",
                "Q_J",
                "observed_order",
                "L2",
                "Linf",
            ]
        )

        def formatted(v):
            return (
                "unavailable"
                if v is None
                else f"{v:.5g}"
                if isinstance(v, (float, np.floating))
                else str(v)
            )

        print(
            f"\n{name.capitalize()} (duration={config['end_time']:g} normalised; viscosity={config['viscosity']:g} m²/s)"
        )
        if name == "convergence":
            print(
                f"Local field errors in {rows[0]['field_units']}; reference: each formulation's finest grid"
            )
        else:
            print(
                "Pressure corrected residual: unavailable (legacy polarisation/full Phi boundary work incomplete)"
            )
        print(" | ".join(fields))
        for row in rows:
            print(" | ".join(formatted(row[k]) for k in fields))
        if latex:
            labels = {
                "final_drift_J": r"$\Delta E$ [J]",
                "max_relative_drift": r"$\max|\Delta E|/E_0$",
                "final_corrected_residual_J": "$R$ [J]",
                "Q_J": r"$K_\perp$ [J]",
                "observed_order": "$p$",
                "L2": "$L_2$ [" + rows[0].get("field_units", "") + "]",
                "Linf": r"$L_\infty$ [" + rows[0].get("field_units", "") + "]",
            }

            def escape(value):
                return formatted(value).replace("_", r"\_")

            text = "\\begin{tabular}{" + "l" * len(fields) + "}\n"
            text += " & ".join(labels.get(k, k) for k in fields) + " \\\\\n\\hline\n"
            for row in rows:
                text += " & ".join(escape(row[k]) for k in fields) + " \\\\\n"
            text += "\\end{tabular}\n"
            (study / (name + "-table.tex")).write_text(text)


def energy_figures(study, config, boundary, groups, *, energy_format="png"):
    """Draw the same energy histories in a combined image or standalone panels."""
    import matplotlib.pyplot as plt

    panels = (
        ("ion-thermal", "Ion thermal change [J]"),
        ("electron-thermal", "Electron thermal change [J]"),
        ("perpendicular-kinetic", "Perpendicular kinetic change [J]"),
        ("shared-relative-change", "Shared FV energy: relative change"),
        (
            "model-relative-change",
            "Model energy: relative change\n(thin: boundary balance residual)",
        ),
        ("estimator-difference", "Shared FV energy − model energy [J]"),
    )
    separate = energy_format == "pdf-panels"
    if separate:
        canvases = [plt.subplots(figsize=(6, 4.5)) for _ in panels]
        axes = [axis for _, axis in canvases]
    else:
        fig, grid = plt.subplots(3, 2, figsize=(11, 10))
        canvases = [(fig, grid)]
        axes = list(grid.flat)

    colours = {"pressure": "tab:blue", "total_energy": "tab:orange"}
    styles = ["--", "-.", ":"]
    for form, rows in groups.items():
        for level, (key, res) in enumerate(rows):
            inv, rep = res["inventory"], res["report"]
            t = np.asarray(inv.t) * 1e3
            kw = {
                "color": colours[form],
                "linestyle": "-"
                if level == len(rows) - 1
                else styles[(len(rows) - 2 - level) % len(styles)],
                "alpha": 0.5 + 0.5 * (level + 1) / len(rows),
                "label": f"{form.replace('_', ' ')} {key[2]}×{key[3]}",
            }
            for axis, name in zip(
                axes[:3], ("ion_thermal", "electron_thermal", "ion_perpendicular")
            ):
                axis.plot(t, inv[name] - inv[name].isel(t=0), **kw)
            axes[3].plot(
                t,
                (inv.plasma_reconstructed - inv.plasma_reconstructed.isel(t=0))
                / inv.plasma_reconstructed.isel(t=0),
                **kw,
            )
            axes[4].plot(
                t,
                (inv.plasma_native - inv.plasma_native.isel(t=0))
                / inv.plasma_native.isel(t=0),
                **kw,
            )
            axes[5].plot(t, inv.plasma_reconstructed - inv.plasma_native, **kw)
            if rep["corrected_residual_J"] is not None:
                axes[4].plot(
                    t,
                    np.asarray(rep["corrected_residual_J"])
                    / float(inv.plasma_native.isel(t=0)),
                    color=colours[form],
                    linestyle=kw["linestyle"],
                    linewidth=0.8,
                    alpha=0.5,
                )
    for axis, (_, title) in zip(axes, panels):
        axis.set_title(title)
        axis.set_xlabel("Time [ms]")
        if separate:
            axis.legend(fontsize=8, ncol=2)

    context = f"{boundary} phi; viscosity {config['viscosity']:g} m²/s"
    residual_note = (
        "Corrected pressure residual unavailable: "
        "incomplete legacy boundary/exchange coverage"
    )
    if separate:
        for (fig, _), (name, _) in zip(canvases, panels):
            title = f"{context}\nparallel kinetic energy = 0"
            if name == "model-relative-change":
                title += "\n" + residual_note.replace(": ", ":\n")
            fig.suptitle(title)
            fig.savefig(study / f"energy-{boundary}-{name}.pdf")
            plt.close(fig)
    else:
        axes[0].legend(fontsize=6, ncol=2)
        fig.suptitle(f"{context}; parallel kinetic energy = 0\n{residual_note}")
        fig.savefig(study / f"energy-{boundary}.png", dpi=160)
        plt.close(fig)


def figures(study, config, results, style=None, *, energy_format="png"):
    os.environ.setdefault("MPLCONFIGDIR", str(study / ".matplotlib"))
    import matplotlib

    matplotlib.use("pdf")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    rc = {
        "font.size": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.grid": False,
        "figure.constrained_layout.use": True,
    }
    saved_style = study / "figure-config.json"
    if saved_style.exists():
        rc.update(json.loads(saved_style.read_text())["rcParams"])
    if style:
        rc.update(json.loads(style.read_text()))
    plt.rcParams.update(rc)
    save_json(
        study / "figure-config.json",
        {
            "rcParams": rc,
            "field": config["field"],
            "snapshot_times_normalised": config["snapshot_times"],
            "convergence_time_normalised": config["convergence_time"],
            "energy_format": energy_format,
        },
    )
    colours = {"pressure": "tab:blue", "total_energy": "tab:orange"}
    for boundary in sorted({r[0] for r in results}):
        groups = {
            f: sorted(
                [
                    (key, res)
                    for key, res in results.items()
                    if key[0] == boundary and key[1] == f
                ],
                key=lambda x: x[0][2],
            )
            for f in colours
        }
        finest = {f: rows[-1][1] for f, rows in groups.items()}
        internal_energy = config["field"] == "ion_pressure"
        upper_arrays = finest["pressure"]["snapshots"]
        if internal_energy:
            # Pressure in Pa and energy density in J/m^3 have the same SI scale.
            upper_arrays = [1.5 * a for a in upper_arrays]
        state_rows = (
            (
                finest["pressure"],
                upper_arrays,
                "Internal energy" if internal_energy else "Pressure",
                r"$3P_i/2$ [J m$^{-3}$]"
                if internal_energy
                else f"{config['field'].replace('_', ' ')} [{finest['pressure']['unit']}]",
            ),
            (
                finest["total_energy"],
                finest["total_energy"]["kinetic_snapshots"],
                r"$K_\perp$",
                r"$K_\perp$ [J m$^{-3}$]",
            ),
        )
        fig, ax = plt.subplots(2, 4, figsize=(12, 5), sharex=True, sharey=True)
        for row, (res, arrays, row_label, colour_label) in enumerate(state_rows):
            vmin, vmax = min(a.min() for a in arrays), max(a.max() for a in arrays)
            signed = row == 0 and config["field"] == "vorticity"
            if signed:
                vmax = max(abs(vmin), abs(vmax))
                vmin = -vmax
            norm = Normalize(vmin=vmin, vmax=vmax)
            for col, a in enumerate(arrays):
                im = ax[row, col].pcolormesh(
                    res["x"] * 1e2,
                    res["z"] * 1e2,
                    a.T,
                    norm=norm,
                    cmap="RdBu_r" if signed else "viridis",
                    rasterized=True,
                    shading="nearest",
                )
                if row == 0:
                    ax[row, col].set_title(
                        rf"$t = {config['snapshot_times'][col]:g}/\omega_{{ci}}$"
                    )
                if row == 1:
                    ax[row, col].set_xlabel("X [cm]")
            ax[row, 0].set_ylabel(f"{row_label}\nZ [cm]")
            fig.colorbar(im, ax=ax[row, :], label=colour_label)
        fig.savefig(study / f"state-{boundary}.pdf")
        fig.savefig(study / f"state-{boundary}.png", dpi=160)
        plt.close(fig)

        energy_figures(study, config, boundary, groups, energy_format=energy_format)

        fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
        for form, rows in groups.items():
            h = np.asarray([res["lx"] / key[2] for key, res in rows])
            q = np.asarray([res["Q"] for key, res in rows])
            seq = global_differences(q, config["refinement"])
            kw = {
                "color": colours[form],
                "marker": "o",
                "label": form.replace("_", " "),
            }
            axes[0].plot(h * 1e3, q, **kw)
            differences = abs(q[:-1] - q[1:])
            mask = differences > 0
            axes[1].loglog(h[:-1][mask] * 1e3, differences[mask], **kw)
            valid = [i for i, v in enumerate(seq) if v["observed_order"] is not None]
            if valid:
                axes[2].plot(
                    h[valid] * 1e3, [seq[i]["observed_order"] for i in valid], **kw
                )
        for axis in axes:
            axis.set_xlabel("Physical ΔX [mm]")
        axes[0].set_ylabel("Common K⊥ [J]")
        axes[1].set_ylabel("Adjacent |ΔK⊥| [J]")
        axes[2].set_ylabel("Observed order")
        axes[0].legend(fontsize=8)
        fig.suptitle(
            f"{boundary} phi; t={config['convergence_time']:g} normalised; orders omitted for oscillatory/roundoff sequences"
        )
        fig.savefig(study / f"global-convergence-{boundary}.pdf")
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
        for form, rows in groups.items():
            for axis, name in zip(axes, ("L2", "Linf")):
                values = [
                    (res["lx"] / key[2] * 1e3, res["error"][name])
                    for key, res in rows[:-1]
                    if res["error"][name] > 0
                ]
                if values:
                    h, e = np.asarray(values).T
                    axis.loglog(
                        h, e, "o-", color=colours[form], label=form.replace("_", " ")
                    )
                axis.set_xlabel("Physical ΔX [mm]")
                axis.set_ylabel(f"{name} [{rows[0][1]['unit']}]")
        axes[0].legend(fontsize=8)
        fig.suptitle(
            f"{boundary} phi; {config['field'].replace('_', ' ')} at t={config['convergence_time']:g}\nEach formulation uses its own finest reference; all physical cells"
        )
        fig.savefig(study / f"local-error-{boundary}.pdf")
        plt.close(fig)


def analyse(study, manifest, *, latex=False, style=None, energy_format="png"):
    config = manifest["configuration"]
    results, energy_rows, conv_rows = {}, [], []
    for row in manifest["runs"]:
        directory = study / row["directory"]
        if digest(directory / "BOUT.inp") != row["input_sha256"]:
            raise AuditError(f"Input was modified: {directory}")
        print(f"Analysing {row['directory']}", flush=True)
        res = analyse_run(directory, config)
        results[(row["boundary"], row["formulation"], row["nx"], row["nz"])] = res
        scalar = {
            k: v for k, v in res["report"].items() if not isinstance(v, (dict, list))
        }
        energy_rows.append(
            {
                "boundary": row["boundary"],
                "formulation": row["formulation"],
                "grid": f"{row['nx']}x{row['nz']}",
                **scalar,
            }
        )
    volumes = [res["report"]["volume_m3"] for res in results.values()]
    if not np.allclose(volumes, volumes[0], rtol=1e-10):
        raise AuditError("Represented physical volumes differ")
    for boundary, form in sorted({key[:2] for key in results}):
        rows = sorted(
            [(key, res) for key, res in results.items() if key[:2] == (boundary, form)],
            key=lambda item: item[0][2],
        )
        fine = rows[-1][1]
        seq = global_differences([res["Q"] for _, res in rows], config["refinement"])
        for i, (key, res) in enumerate(rows):
            errors = {k: None for k in ("L2", "Linf", "relative_L2", "relative_Linf")}
            if i < len(rows) - 1:
                reference = restrict(
                    fine["comparison"],
                    fine["x"],
                    fine["z"],
                    res["x"],
                    res["z"],
                    lx=res["lx"],
                    lz=res["lz"],
                    representation="point"
                    if config["field"] == "vorticity"
                    else "cell_average",
                )
                errors = error_norms(res["comparison"], reference, res["volume"])
            res["error"] = errors
            conv_rows.append(
                {
                    "boundary": boundary,
                    "formulation": form,
                    "grid": f"{key[2]}x{key[3]}",
                    "spacing_x_m": res["lx"] / key[2],
                    "spacing_z_m": res["lz"] / key[3],
                    "field_units": res["unit"],
                    **seq[i],
                    **errors,
                }
            )
    write_tables(study, energy_rows, conv_rows, config, latex=latex)
    figures(study, config, results, style, energy_format=energy_format)
    return results


def pilot(study, manifest, exe, timeout, results):
    """Measure temporal/solver/cadence sensitivity against spatial differences."""
    config = manifest["configuration"]
    evidence = []
    for row in manifest["runs"]:
        if runtime_artifacts(exe) != manifest["runtime_artifacts"]:
            raise AuditError(
                "Executable or model libraries changed during the accuracy pilot"
            )
        if row["level"] < config["levels"] - 2:
            continue
        key = (row["boundary"], row["formulation"], row["nx"], row["nz"])
        base = results[key]
        for variant in ("tight", "cadence"):
            c = dict(config)
            if variant == "tight":
                for name in ("atol", "rtol", "phi_atol", "phi_rtol", "vort_atol"):
                    c[name] /= 10
                c["max_timestep"] /= 2
            else:
                c["outputs"] *= 2
            directory = study / "pilot" / f"{row['directory']}-{variant}"
            text = input_text(c, row)
            if not (directory / "status.json").exists():
                directory.mkdir(parents=True, exist_ok=True)
                (directory / "BOUT.inp").write_text(text)
                build_startup(exe, directory, c, timeout)
                execute(exe, directory, c, timeout)
            elif (directory / "BOUT.inp").read_text() != text:
                raise AuditError("Pilot input hash/configuration mismatch")
            res = analyse_run(directory, c)
            q_error = abs(res["Q"] - base["Q"])
            group = sorted(
                [(k, r) for k, r in results.items() if k[:2] == key[:2]],
                key=lambda x: x[0][2],
            )
            spatial = min(
                abs(group[i][1]["Q"] - group[i + 1][1]["Q"])
                for i in range(len(group) - 1)
            )
            field_error = error_norms(
                res["comparison"], base["comparison"], base["volume"]
            )["L2"]
            errors = [
                r["error"]["L2"] for _, r in group[:-1] if r["error"]["L2"] is not None
            ]
            residual_error = None
            if base["report"]["corrected_residual_J"] is not None:
                a = (
                    np.asarray(res["report"]["corrected_residual_J"])[::2]
                    if variant == "cadence"
                    else np.asarray(res["report"]["corrected_residual_J"])
                )
                residual_error = float(
                    abs(a - base["report"]["corrected_residual_J"]).max()
                )
            evidence.append(
                {
                    "directory": str(directory.relative_to(study)),
                    "variant": variant,
                    "Q_difference_J": q_error,
                    "minimum_spatial_Q_difference_J": spatial,
                    "Q_ratio_to_spatial_difference": q_error / spatial
                    if spatial
                    else None,
                    "L2_temporal_difference": field_error,
                    "L2_ratio_to_spatial_error": field_error / min(errors)
                    if errors and min(errors) > 0
                    else None,
                    "max_residual_difference_J": residual_error,
                    "residual_ratio_to_spatial_Q_difference": residual_error / spatial
                    if residual_error is not None and spatial
                    else None,
                    "estimated_fine_power_sampling_error_J": residual_error / 3
                    if variant == "cadence" and residual_error is not None
                    else None,
                    "target_under_10_percent": bool(
                        spatial > 0
                        and q_error < 0.1 * spatial
                        and errors
                        and field_error < 0.1 * min(errors)
                        and (residual_error is None or residual_error < 0.1 * spatial)
                    ),
                }
            )
    save_json(study / "accuracy-pilot.json", evidence)
    return evidence


def grid_rows(config):
    boundaries = (
        ("fixed", "relax")
        if config["phi_boundary"] == "both"
        else (config["phi_boundary"],)
    )
    for boundary in boundaries:
        for formulation in ("pressure", "total_energy"):
            for level in range(config["levels"]):
                nx = config["base_nx"] * config["refinement"] ** level
                nz = config["base_nz"] * config["refinement"] ** level
                yield {
                    "boundary": boundary,
                    "formulation": formulation,
                    "level": level,
                    "nx": nx,
                    "nz": nz,
                    "directory": f"{boundary}-{formulation}-{nx}x{nz}",
                }


def add_inputs(study, config, rows):
    # Preflight all collisions before creating anything. A previous interrupted
    # extension can leave matching inputs that are safe to adopt.
    for row in rows:
        path = study / row["directory"] / "BOUT.inp"
        if path.exists() and path.read_text() != input_text(config, row):
            raise AuditError(f"Existing input differs: {path}")
    for row in rows:
        directory = study / row["directory"]
        directory.mkdir(exist_ok=True)
        path = directory / "BOUT.inp"
        if not path.exists():
            path.write_text(input_text(config, row))
        row["input_sha256"] = digest(path)


def compatible_configuration(saved, config):
    return config["levels"] >= saved["levels"] and {
        k: v for k, v in saved.items() if k != "levels"
    } == {k: v for k, v in config.items() if k != "levels"}


def extend_study(study, manifest, config):
    old = manifest["configuration"]
    if not compatible_configuration(old, config):
        raise AuditError(
            "Extension may only increase --levels; other settings must match"
        )
    if config["levels"] == old["levels"]:
        return
    for row in manifest["runs"]:
        if digest(study / row["directory"] / "BOUT.inp") != row["input_sha256"]:
            raise AuditError(f"Modified input: {row['directory']}")
    rows = [row for row in grid_rows(config) if row["level"] >= old["levels"]]
    add_inputs(study, config, rows)
    manifest.setdefault("extensions", []).append(
        {
            "previous_levels": old["levels"],
            "levels": config["levels"],
            "added_cases": [row["directory"] for row in rows],
            "recorded_utc": datetime.now(timezone.utc).isoformat(),
        }
    )
    manifest["configuration"] = config
    manifest["runs"].extend(rows)
    save_json(study / "manifest.json", manifest)
    print(
        f"Extended study from {old['levels']} to {config['levels']} levels", flush=True
    )


def select_cases(manifest, names):
    rows = manifest["runs"]
    if not names:
        return rows
    unknown = set(names) - {row["directory"] for row in rows}
    if unknown:
        raise AuditError(f"Unknown case directories: {', '.join(sorted(unknown))}")
    return [row for row in rows if row["directory"] in names]


def check_prepared_files(directory, case):
    if case["schema_version"] != 1:
        raise AuditError("Unsupported remote case schema")
    for name, expected in case["input_files"].items():
        if Path(name).name != name or digest(directory / name) != expected:
            raise AuditError(f"Prepared input changed: {directory / name}")


def archive_attempt(directory):
    archive = directory / f"failed-{time.time_ns()}"
    archive.mkdir()
    files = list(directory.glob("BOUT.*.nc")) + list(directory.glob("BOUT.log*"))
    files += [
        directory / name
        for name in (
            "initial-state.nc",
            "startup",
            "status.json",
            "run.log",
            "remote-case.json",
            "remote-execution.json",
            "run-remote.py",
        )
    ]
    for file in files:
        if file.exists():
            file.rename(archive / file.name)


def ensure_startup(exe, directory, config, timeout):
    prepared = directory / "remote-case.json"
    if prepared.exists() and not (directory / "status.json").exists():
        if (
            list(directory.glob("BOUT.dmp.*.nc"))
            or (directory / "remote-execution.json").exists()
        ):
            raise AuditError(
                f"Unregistered outputs in {directory}; use --import-runs first"
            )
        check_prepared_files(directory, json.loads(prepared.read_text()))
        return
    if (directory / "startup").exists():
        archive_attempt(directory)
    if not (directory / "initial-state.nc").exists() or not list(
        directory.glob("BOUT.restart.*.nc")
    ):
        build_startup(exe, directory, config, timeout)


def prepare_remote(directory, row, manifest):
    files = [directory / "BOUT.inp", directory / "initial-state.nc"]
    files += sorted(directory.glob("BOUT.restart.*.nc"))
    if len(files) != manifest["configuration"]["nproc"] + 2:
        raise AuditError(f"Missing matched restart files: {directory}")
    path = directory / "remote-case.json"
    if not path.exists():
        save_json(
            path,
            {
                "schema_version": 1,
                "study": manifest["created_utc"],
                "case": row,
                "configuration": manifest["configuration"],
                "input_files": {file.name: digest(file) for file in files},
                "expected_revisions": manifest["revisions"],
                "prepared_utc": datetime.now(timezone.utc).isoformat(),
            },
        )
    else:
        check_prepared_files(directory, json.loads(path.read_text()))
    shutil.copyfile(HERE / "publication_remote.py", directory / "run-remote.py")
    print(f"Prepared {row['directory']} for local or remote execution", flush=True)


def link_or_copy(source, target):
    try:
        os.link(source, target)
    except OSError:
        shutil.copy2(source, target)


def validate_supplied_options(directory, options):
    from boutdata.data import BoutOptionsFile

    expected = BoutOptionsFile(str(directory / "BOUT.inp")).as_dict()

    def compare(given, actual, prefix=""):
        for name, value in given.items():
            key = f"{prefix}:{name}" if prefix else name
            if name not in actual:
                raise AuditError(f"Missing supplied option in returned settings: {key}")
            if isinstance(value, dict):
                if not isinstance(actual[name], dict):
                    raise AuditError(f"Returned option section differs: {key}")
                compare(value, actual[name], key)
            elif actual[name] != value:
                raise AuditError(f"Returned setting differs from supplied input: {key}")

    compare(expected, options.as_dict())


def validate_import(directory, row, config):
    check = verify_completed_outputs(directory, config)
    run = open_run(directory)
    try:
        validate_supplied_options(directory, run.data.options)
        if (
            run.formulation != row["formulation"]
            or run.data.sizes["x"] != row["nx"] + 4
            or run.data.sizes["z"] != row["nz"]
            or run.data.sizes["y"] != 1
        ):
            raise AuditError(
                f"Imported grid/formulation does not match {row['directory']}"
            )
        coordinates(run)
        validate_times(run, config)
        initial = validate_initial_state(directory, run)
    finally:
        run.data.close()
    return check, initial


def import_execution(source, destination, row, manifest):
    path = source / "remote-execution.json"
    if not path.exists():
        print(
            f"Warning: no remote-execution.json for {row['directory']}; "
            "remote exit code and build provenance are unknown.",
            file=sys.stderr,
        )
        return None
    execution = json.loads(path.read_text())
    prepared = destination / "remote-case.json"
    if (
        execution.get("schema_version") != 1
        or execution.get("case") != row["directory"]
        or not prepared.exists()
        or execution.get("prepared_case_sha256") != digest(prepared)
        or execution.get("input_files")
        != json.loads(prepared.read_text())["input_files"]
    ):
        raise AuditError(f"Remote execution does not match prepared case: {source}")
    if execution.get("runtime_artifacts") != execution.get("runtime_artifacts_after"):
        raise AuditError(
            f"Remote executable/libraries changed during execution: {source}"
        )
    for name in ("hermes", "bout"):
        expected = manifest.get("revisions", {}).get(name, {}).get("revision")
        actual = execution.get("revisions", {}).get(name, {}).get("revision")
        if expected and actual and expected != actual:
            raise AuditError(f"Remote {name} revision differs from study: {source}")
    return execution


def import_runs(study, manifest, names=None, source_root=None):
    """Stage, validate and register returned cases without requiring Hermes."""
    rows = select_cases(manifest, names)
    source_root = source_root.resolve() if source_root else study
    imported = 0
    for row in rows:
        destination = study / row["directory"]
        source = source_root / row["directory"]
        if completed_run(destination, manifest["configuration"]):
            print(f"Reusing complete {row['directory']}", flush=True)
            continue
        files = sorted(source.glob("BOUT.dmp.*.nc"))
        if not files and not names:
            continue
        if (
            digest(source / "BOUT.inp") != row["input_sha256"]
            or digest(destination / "BOUT.inp") != row["input_sha256"]
        ):
            raise AuditError(f"Modified returned/local input: {source}")
        execution = import_execution(source, destination, row, manifest)
        with tempfile.TemporaryDirectory(prefix=".import-", dir=study) as staging:
            stage = Path(staging)
            for file in files:
                link_or_copy(file, stage / file.name)
            for name in ("BOUT.inp", "initial-state.nc"):
                link_or_copy(destination / name, stage / name)
            if (source / "BOUT.settings").is_file():
                link_or_copy(source / "BOUT.settings", stage / "BOUT.settings")
            check, initial = validate_import(stage, row, manifest["configuration"])
            status = {
                "origin": "external",
                "exit_code": execution["exit_code"] if execution else None,
                "outputs_complete": True,
                "output_check": check,
                "initial_state_check": initial,
                "imported_utc": datetime.now(timezone.utc).isoformat(),
                "source": str(source),
                "remote_execution": execution,
                "options_source": "BOUT.settings"
                if (stage / "BOUT.settings").exists()
                else "BOUT.inp",
                "files": {file.name: file.stat().st_size for file in files},
            }
            old_status = destination / "status.json"
            if old_status.exists():
                status["previous_status"] = json.loads(old_status.read_text())
            if source != destination:
                old_files = list(destination.glob("BOUT.dmp.*.nc"))
                old_files += list(destination.glob("BOUT.log*"))
                old_files += [
                    destination / name
                    for name in (
                        "run.log",
                        "BOUT.settings",
                        "remote-execution.json",
                        "status.json",
                    )
                    if (destination / name).exists()
                ]
                if old_files:
                    archive = destination / f"pre-import-{time.time_ns()}"
                    archive.mkdir()
                    for file in old_files:
                        file.rename(archive / file.name)
                for file in files:
                    (stage / file.name).replace(destination / file.name)
                for pattern in (
                    "BOUT.log*",
                    "run.log",
                    "remote-execution.json",
                    "BOUT.settings",
                ):
                    for file in source.glob(pattern):
                        if file.is_file():
                            shutil.copy2(file, destination / file.name)
            save_json(old_status, status)
        manifest.setdefault("imports", []).append(
            {
                "case": row["directory"],
                "imported_utc": status["imported_utc"],
                "source": str(source),
                "remote_execution": execution,
                "runner_sha256": digest(__file__),
            }
        )
        save_json(study / "manifest.json", manifest)
        imported += 1
        print(f"Imported complete {row['directory']}", flush=True)
    print(
        f"Imported {imported} cases. Resume to run remaining cases or use --analyse-only."
    )
    return imported


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if args.extend:
        saved = json.loads((args.extend / "manifest.json").read_text())
        p.set_defaults(**saved["configuration"])
        if "path" in saved["executable"]:
            p.set_defaults(exe=Path(saved["executable"]["path"]))
        args = p.parse_args(argv)
        if args.levels <= saved["configuration"]["levels"]:
            raise AuditError(
                "--extend requires --levels greater than the saved level count"
            )
        args.resume = args.extend
    if args.import_from and not args.import_runs:
        raise AuditError("--import-from requires --import-runs")
    if args.prepare_only and (args.analyse_only or args.import_runs):
        raise AuditError("--prepare-only requires a new, resumed or extended study")
    if args.cases and args.analyse_only:
        raise AuditError("--cases cannot restrict --analyse-only")
    if args.allow_runner_change and not args.resume:
        raise AuditError("--allow-runner-change requires --resume")
    if args.import_runs:
        study = args.import_runs.resolve()
        manifest = json.loads((study / "manifest.json").read_text())
        if manifest["schema_version"] != SCHEMA:
            raise AuditError("Unsupported manifest schema")
        import_runs(study, manifest, args.cases, args.import_from)
        return 0
    if args.analyse_only:
        study = args.analyse_only.resolve()
        manifest = json.loads((study / "manifest.json").read_text())
        if manifest["schema_version"] != SCHEMA:
            raise AuditError("Unsupported manifest schema")
        analyse(
            study,
            manifest,
            latex=args.latex,
            style=args.style,
            energy_format=args.energy_format,
        )
        print(f"Figures and tables: {study}")
        return 0
    config = configuration(args)
    exe = args.exe.resolve()
    if not exe.is_file():
        raise AuditError(f"Executable not found: {exe}")
    if args.timeout <= 0:
        raise AuditError("Timeout must be positive")
    exe_hash = digest(exe)
    runtime_hashes = runtime_artifacts(exe)
    template_hash = digest(HERE / "publication-input.inp.template")
    if args.resume:
        study = args.resume.resolve()
        manifest = json.loads((study / "manifest.json").read_text())
        if (
            not compatible_configuration(manifest["configuration"], config)
            or manifest["executable"]["sha256"] != exe_hash
            or manifest["template_sha256"] != template_hash
            or manifest.get("runtime_artifacts") != runtime_hashes
        ):
            raise AuditError(
                "Resume requires matching configuration, executable and template hashes"
            )
        current_runner = digest(__file__)
        previous_runner = manifest.get(
            "resume_runner_sha256", manifest["runner_sha256"]
        )
        if previous_runner != current_runner:
            if not args.allow_runner_change:
                raise AuditError(
                    "Runner changed; use --allow-runner-change with --resume to record "
                    "and accept the updated runner"
                )
            manifest.setdefault("runner_updates", []).append(
                {
                    "previous_sha256": previous_runner,
                    "sha256": current_runner,
                    "recorded_utc": datetime.now(timezone.utc).isoformat(),
                }
            )
            manifest["resume_runner_sha256"] = current_runner
            save_json(study / "manifest.json", manifest)
        extend_study(study, manifest, config)
    else:
        args.output_root.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        study = args.output_root.resolve() / f"interchange-{stamp}"
        study.mkdir()
        from importlib.metadata import version

        rows = list(grid_rows(config))
        add_inputs(study, config, rows)
        manifest = {
            "schema_version": SCHEMA,
            "created_utc": stamp,
            "configuration": config,
            "runs": rows,
            "executable": {"path": str(exe), "sha256": exe_hash},
            "runtime_artifacts": runtime_hashes,
            "template_sha256": template_hash,
            "runner_sha256": digest(__file__),
            "formulation_differences": DIFFERENCES,
            "revisions": {
                "hermes": revision(SOURCE_ROOT),
                "bout": revision(SOURCE_ROOT / "external/BOUT-dev"),
                "energy_analysis": revision(ROOT),
            },
            "environment": {
                "python": sys.executable,
                "version": platform.python_version(),
                "numpy": np.__version__,
                "xarray": xr.__version__,
                "scipy": version("scipy"),
                "matplotlib": version("matplotlib"),
                "netCDF4": version("netCDF4"),
                "MPI_launcher": shutil.which("mpiexec"),
                "MPI_nproc": config["nproc"],
            },
            "field_representation": "cell averages for pressure/density; cubic-X/Fourier-Z samples for vorticity",
        }
        save_json(study / "manifest.json", manifest)
    print(f"Study: {study}", flush=True)
    for row in select_cases(manifest, args.cases):
        directory = study / row["directory"]
        if runtime_artifacts(exe) != runtime_hashes:
            raise AuditError(
                "Executable or model libraries changed during the simulation matrix"
            )
        if digest(directory / "BOUT.inp") != row["input_sha256"]:
            raise AuditError(f"Modified input: {directory}")
        if completed_run(directory, config):
            print(f"Reusing complete {row['directory']}", flush=True)
            continue
        ensure_startup(exe, directory, config, args.timeout)
        if args.prepare_only:
            prepare_remote(directory, row, manifest)
            continue
        print(f"Running {row['directory']}", flush=True)
        execute(exe, directory, config, args.timeout)
    if args.prepare_only:
        print(f"Prepared cases: {study}. No full simulations or analysis were run.")
        return 0
    if args.cases and any(
        not completed_run(study / row["directory"], config) for row in manifest["runs"]
    ):
        print("Selected cases finished. Resume to complete the matrix before analysis.")
        return 0
    results = analyse(
        study,
        manifest,
        latex=args.latex,
        style=args.style,
        energy_format=args.energy_format,
    )
    if config["accuracy_pilot"]:
        evidence = pilot(study, manifest, exe, args.timeout, results)
        print(
            f"Accuracy pilot: {sum(v['target_under_10_percent'] for v in evidence)}/{len(evidence)} comparisons below 10% spatial error",
            flush=True,
        )
    print(f"Figures and tables: {study}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print(
            "Publication study interrupted; outputs and status are retained for resume.",
            file=sys.stderr,
        )
        sys.exit(130)
    except (AuditError, FileNotFoundError) as error:
        print(f"Publication study stopped: {error}", file=sys.stderr)
        sys.exit(1)
