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
source_path = HERE / "source-root.txt"
SOURCE_ROOT = Path(source_path.read_text().strip()) if source_path.exists() else ROOT
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
from energy_analysis.operators import perpendicular_operator
from energy_analysis.pressure import legacy_operator, legacy_viscous_work
from energy_analysis.report import revision

SCHEMA = 1
DIFFERENCES = [
    "Pressure uses a potential-flow polarisation approximation and electron compression; total energy uses ion-current partition and species transport.",
    "Pressure polarisation potential includes electron pressure with its signed mass/charge weight; total energy includes ions only.",
    "Pressure native inversion uses FD X / spectral Z; total energy uses face/corner FV X-Z.",
    "Pressure retains cell-centred corner-gradient viscous heating; total energy uses conservative face work.",
    "Initial sine pressure has incompatible homogeneous radial Neumann derivatives; boundary-adjacent cells remain in the error norms.",
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
    p.add_argument("--phi-boundary-timescale", type=float, default=1e-4)
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
    p.add_argument("--analyse-only", type=Path, metavar="STUDY")
    p.add_argument("--resume", type=Path, metavar="STUDY")
    p.add_argument("--diagnose", action="store_true")
    p.add_argument("--atol", type=float, default=1e-11)
    p.add_argument("--rtol", type=float, default=1e-8)
    p.add_argument("--phi-rtol", type=float, default=1e-8)
    p.add_argument(
        "--phi-atol", type=float, default=1e-5, help="Potential iteration scale [eV]"
    )
    p.add_argument("--vort-atol", type=float, default=1e-12)
    p.add_argument("--max-timestep", type=float, default=50)
    p.add_argument(
        "--accuracy-pilot",
        action="store_true",
        help="Repeat two finest grids with tightened tolerances/timesteps and output cadence",
    )
    p.add_argument("--latex", action="store_true")
    p.add_argument("--style", type=Path, help="JSON Matplotlib rcParams overrides")
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
            "analyse_only",
            "style",
            "latex",
            "smoke",
            "timeout",
        )
    }
    if c["levels"] < 2 or min(c["base_nx"], c["base_nz"]) < 4 or c["refinement"] < 2:
        raise AuditError("Use >=2 grids, >=4 physical cells, refinement>=2")
    for key in (
        "end_time",
        "phi_boundary_timescale",
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
    c = config
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
        "relax": str(row["boundary"] == "relax").lower(),
        "tau": c["phi_boundary_timescale"],
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
            else f"density_floor = 1e-7\nneumann_boundary_average_z = false\nphi_atol = {c['phi_atol']}\nphi_rtol = {c['phi_rtol']}\nvort_atol = {c['vort_atol']}\nphi_max_iterations = 400"
        ),
        "polarisation_options": (
            "[polarisation_drift]\nboussinesq = true\naverage_atomic_mass = 1\ndiamagnetic_polarisation = true\nadvection = true\ndiagnose = "
            + str(c["diagnose"]).lower()
        )
        if pressure
        else "",
    }
    text = string.Template((HERE / "publication-input.inp").read_text()).substitute(
        values
    )
    # Only the formulation's evolving fields have initial-profile sections.
    import re

    section = "Eion" if pressure else "Ph+"
    text = re.sub(
        r"\[" + re.escape(section) + r"\]\n.*?(?=\n\[)", "", text, flags=re.DOTALL
    )
    if row["boundary"] == "fixed":
        text = re.sub(r"^phi_boundary_timescale = .*\n", "", text, flags=re.MULTILINE)
        text = text.replace("phi_core_averagey = false\n", "")
    return text


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
    save_json(directory / ("startup-status.json" if seed else "status.json"), status)
    if interrupted:
        raise KeyboardInterrupt
    if code:
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
    pressure = 1 + 0.1 * np.sin(2 * np.pi * x[:, None] - z[None, :])
    harmonic = 1 - 0.1 * np.sin(z)[None, :] * np.cosh(
        2 * np.pi * (x[:, None] - 0.5)
    ) / np.cosh(np.pi)
    potential = harmonic - pressure
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
    psi = p * (1 - 1 / 1836) if run.formulation == "pressure" else p
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
    energy = 1.5 * p + op.kperp
    for i in range(2):
        energy[{"x": i}] = energy.isel(x=3 - i)
        energy[{"x": -i - 1}] = energy.isel(x=-4 + i)
    return {
        "Pe": p,
        "Ph+": p,
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
    relaxing = run.data.options[run.model].get_bool("phi_boundary_relax")
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
    expected = np.linspace(0, config["end_time"], config["outputs"] + 1)
    if values.shape != expected.shape or not np.allclose(
        values, expected, rtol=0, atol=config["end_time"] * 1e-10
    ):
        raise AuditError(f"Incomplete/mismatched output times: {run.source}")
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


def analyse_run(directory, config):
    status = json.loads((directory / "status.json").read_text())
    if status["exit_code"] != 0:
        raise AuditError(f"Refusing failed run {directory}")
    run = open_run(directory)
    save_json(directory / "resolved-options.json", run.data.options.as_dict())
    snapshots, convergence = validate_times(run, config)
    with xr.open_dataset(directory / "initial-state.nc") as saved:
        initial = saved.load()
    checks = {}
    first = replace(run, data=run.data.isel(t=[0]).compute())
    for name in ("Pe", "Nh+", "phi", "Vort"):
        checks[name] = float(
            abs(run.interior(first.data[name].isel(t=0) - initial[name])).max()
        )
        if checks[name] > 1e-8:
            raise AuditError(f"Initial {name} mismatch {checks[name]}: {directory}")
    pi = (
        first.data["Ph+"]
        if run.formulation == "pressure"
        else first.data["Nh+"] * first.data.Tion
    )
    checks["ion_pressure"] = float(
        abs(run.interior(pi.isel(t=0) - initial["Ph+"])).max()
    )
    if checks["ion_pressure"] > 1e-8:
        raise AuditError("Initial ion pressure mismatch")
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
    if (
        not run.data.options[run.model].get_bool("phi_boundary_relax")
        and float(abs(bx).max()) > 1e-9
    ):
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


def figures(study, config, results, style=None):
    os.environ.setdefault("MPLCONFIGDIR", str(study / ".matplotlib"))
    import matplotlib

    matplotlib.use("pdf")
    import matplotlib.pyplot as plt

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
        },
    )
    colours = {"pressure": "tab:blue", "total_energy": "tab:orange"}
    styles = ["--", "-.", ":"]
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
        arrays = [a for res in finest.values() for a in res["snapshots"]]
        vmin, vmax = min(a.min() for a in arrays), max(a.max() for a in arrays)
        if config["field"] == "vorticity":
            vmax = max(abs(vmin), abs(vmax))
            vmin = -vmax
        fig, ax = plt.subplots(2, 4, figsize=(12, 5), sharex=True, sharey=True)
        for row, (form, res) in enumerate(finest.items()):
            for col, a in enumerate(res["snapshots"]):
                im = ax[row, col].pcolormesh(
                    res["x"] * 1e2,
                    res["z"] * 1e2,
                    a.T,
                    vmin=vmin,
                    vmax=vmax,
                    cmap="RdBu_r" if config["field"] == "vorticity" else "viridis",
                    rasterized=True,
                    shading="nearest",
                )
                if row == 0:
                    ax[row, col].set_title(f"t = {config['snapshot_times'][col]:g}")
                if row == 1:
                    ax[row, col].set_xlabel("X [cm]")
            ax[row, 0].set_ylabel(f"{form.replace('_', ' ')}\nZ [cm]")
        fig.colorbar(
            im,
            ax=ax,
            label=f"{config['field'].replace('_', ' ')} [{next(iter(finest.values()))['unit']}]",
        )
        fig.suptitle(f"{boundary} phi; finest grid; normalised snapshot times")
        fig.savefig(study / f"state-{boundary}.pdf")
        plt.close(fig)

        fig, axes = plt.subplots(3, 2, figsize=(11, 10))
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
                    axes.flat[:3],
                    ("ion_thermal", "electron_thermal", "ion_perpendicular"),
                ):
                    axis.plot(t, inv[name] - inv[name].isel(t=0), **kw)
                axes.flat[3].plot(
                    t,
                    (inv.plasma_reconstructed - inv.plasma_reconstructed.isel(t=0))
                    / inv.plasma_reconstructed.isel(t=0),
                    **kw,
                )
                axes.flat[4].plot(
                    t,
                    (inv.plasma_native - inv.plasma_native.isel(t=0))
                    / inv.plasma_native.isel(t=0),
                    **kw,
                )
                axes.flat[5].plot(t, inv.plasma_reconstructed - inv.plasma_native, **kw)
                if rep["corrected_residual_J"] is not None:
                    axes.flat[4].plot(
                        t,
                        np.asarray(rep["corrected_residual_J"])
                        / float(inv.plasma_native.isel(t=0)),
                        color=colours[form],
                        linestyle=kw["linestyle"],
                        linewidth=0.8,
                        alpha=0.5,
                    )
        titles = [
            "Ion thermal change [J]",
            "Electron thermal change [J]",
            "Perpendicular kinetic change [J]",
            "Shared FV energy: relative change",
            "Model energy: relative change\n(thin: boundary balance residual)",
            "Shared FV energy − model energy [J]",
        ]
        for axis, title in zip(axes.flat, titles):
            axis.set_title(title)
            axis.set_xlabel("Time [ms]")
        axes.flat[0].legend(fontsize=6, ncol=2)
        fig.suptitle(
            f"{boundary} phi; viscosity {config['viscosity']:g} m²/s; parallel kinetic energy = 0\nCorrected pressure residual unavailable: incomplete legacy boundary/exchange coverage"
        )
        fig.savefig(study / f"energy-{boundary}.pdf")
        plt.close(fig)

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


def analyse(study, manifest, *, latex=False, style=None):
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
    figures(study, config, results, style)
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


def main(argv=None):
    args = parser().parse_args(argv)
    if args.analyse_only and args.resume:
        raise AuditError("Choose analyse-only or resume")
    if args.analyse_only:
        study = args.analyse_only.resolve()
        manifest = json.loads((study / "manifest.json").read_text())
        if manifest["schema_version"] != SCHEMA:
            raise AuditError("Unsupported manifest schema")
        analyse(study, manifest, latex=args.latex, style=args.style)
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
    template_hash = digest(HERE / "publication-input.inp")
    if args.resume:
        study = args.resume.resolve()
        manifest = json.loads((study / "manifest.json").read_text())
        if (
            manifest["configuration"] != config
            or manifest["executable"]["sha256"] != exe_hash
            or manifest["template_sha256"] != template_hash
            or manifest["runner_sha256"] != digest(__file__)
            or manifest.get("runtime_artifacts") != runtime_hashes
        ):
            raise AuditError(
                "Resume requires matching configuration, executable, template and runner hashes"
            )
    else:
        args.output_root.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        study = args.output_root.resolve() / f"interchange-{stamp}"
        study.mkdir()
        from importlib.metadata import version

        import xhermes

        boundaries = (
            ("fixed", "relax")
            if config["phi_boundary"] == "both"
            else (config["phi_boundary"],)
        )
        rows = []
        for b in boundaries:
            for f in ("pressure", "total_energy"):
                for level in range(config["levels"]):
                    nx = config["base_nx"] * config["refinement"] ** level
                    nz = config["base_nz"] * config["refinement"] ** level
                    row = {
                        "boundary": b,
                        "formulation": f,
                        "level": level,
                        "nx": nx,
                        "nz": nz,
                        "directory": f"{b}-{f}-{nx}x{nz}",
                    }
                    directory = study / row["directory"]
                    directory.mkdir()
                    (directory / "BOUT.inp").write_text(input_text(config, row))
                    row["input_sha256"] = digest(directory / "BOUT.inp")
                    rows.append(row)
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
                "xhermes": revision(Path(xhermes.__file__).resolve().parents[1]),
            },
            "environment": {
                "python": sys.executable,
                "version": platform.python_version(),
                "xhermes": version("xhermes"),
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
    for row in manifest["runs"]:
        directory = study / row["directory"]
        if runtime_artifacts(exe) != runtime_hashes:
            raise AuditError(
                "Executable or model libraries changed during the simulation matrix"
            )
        if digest(directory / "BOUT.inp") != row["input_sha256"]:
            raise AuditError(f"Modified input: {directory}")
        status = directory / "status.json"
        if status.exists() and json.loads(status.read_text())["exit_code"] == 0:
            print(f"Reusing successful {row['directory']}", flush=True)
            continue
        if (directory / "startup").exists() and (
            not status.exists() or json.loads(status.read_text())["exit_code"] != 0
        ):
            # Archive failed attempts; never pool their dumps with this attempt.
            archive = directory / f"failed-{time.time_ns()}"
            archive.mkdir()
            for file in directory.glob("BOUT.*.nc"):
                file.rename(archive / file.name)
            if (directory / "initial-state.nc").exists():
                (directory / "initial-state.nc").rename(archive / "initial-state.nc")
            if (directory / "startup").exists():
                (directory / "startup").rename(archive / "startup")
        print(f"Running {row['directory']}", flush=True)
        if not (directory / "initial-state.nc").exists() or not list(
            directory.glob("BOUT.restart.*.nc")
        ):
            build_startup(exe, directory, config, args.timeout)
        execute(exe, directory, config, args.timeout)
    results = analyse(study, manifest, latex=args.latex, style=args.style)
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
