"""Bounded COMPASS audits; the general tokamak power audit remains incomplete.

Read each rank with both processor halos intact. X-Z operators have no Y
derivatives, and parallel reconstruction uses those halos without flattening
the core/private-flux connections into an ordinary periodic array.
"""

import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import xarray as xr
import xbout  # noqa: F401  Registers dataset.options.
from boutdata.data import BoutOptionsFile
from ion_energy_balance import ROOT, execute, histories, require, write_json

sys.path.insert(0, str(ROOT / "tools"))
from energy_analysis import energy_balance
from energy_analysis.dataset import EnergyRun, Species
from energy_analysis.inventory import (
    compare_diagnostics,
    energy_densities,
    energy_inventory,
    soft_floor,
)
from energy_analysis.operators import perpendicular_operator
from energy_analysis.transport import parallel_flux

MODEL = "evolve_vorticity_ion_energy"
ESSENTIAL = (
    "phi",
    "Vort",
    "Eion",
    "Tion",
    "Pe",
    "Ne",
    "Nd+",
    "Nimpurity",
    "NVd+",
    "NVimpurity",
)


def read_ranks(folder):
    options = BoutOptionsFile(
        str(folder / "BOUT.settings"), gridfilename=str(folder / "compass-36x48.grd.nc")
    )
    runs = []
    for path in sorted(folder.glob("BOUT.dmp.*.nc")):
        with xr.open_dataset(path) as saved:
            ds = saved.load()
        meta = {
            n: ds[n].item()
            for n in (
                "MXG",
                "MYG",
                "PE_XIND",
                "PE_YIND",
                "NXPE",
                "NYPE",
                "MXSUB",
                "MYSUB",
                "ixseps1",
                "jyseps1_1",
                "jyseps2_2",
                "ny",
                "Nnorm",
                "Tnorm",
                "Omega_ci",
                "rho_s0",
            )
        }
        ds = ds.assign_coords(
            t=ds.t_array.values,
            x=np.arange(ds.sizes["x"]),
            y=np.arange(ds.sizes["y"]),
            z=np.arange(ds.sizes["z"]),
        )
        # BOUT dump values are already normalised; conversion attrs describe SI.
        for name in ds:
            ds[name].attrs["units_type"] = "hermes"
        ds.attrs.update(options=options, metadata=meta)
        species = tuple(
            Species(
                name,
                ds.options[name].evaluate_scalar("AA"),
                ds.options[name].evaluate_scalar("charge"),
                tuple(t.strip() for t in ds.options[name]["type"].split(",")),
            )
            for name in ("e", "d+", "impurity")
        )
        runs.append(
            EnergyRun(
                ds,
                species,
                MODEL,
                int(meta["MXG"]),
                int(meta["MYG"]),
                1.602176634e-19 * meta["Tnorm"] * meta["Nnorm"],
                meta["rho_s0"],
                meta["Omega_ci"],
                str(path),
            )
        )
    return sorted(runs, key=lambda r: r.data.metadata["PE_YIND"])


def boundary_rates(runs):
    """Midpoints and rates, including only the connected core in its Y mean."""
    values = []
    core = []
    meta = runs[0].data.metadata
    for run in runs:
        y = slice(run.myg, -run.myg)
        inner = run.data.phi.isel(x=run.mxg, y=y).mean("z")
        outer = run.data.phi.isel(x=-run.mxg - 1, y=y).mean("z")
        binner = (inner + run.data.phi.isel(x=run.mxg - 1, y=y).mean("z")) / 2
        bouter = (outer + run.data.phi.isel(x=-run.mxg, y=y).mean("z")) / 2
        global_y = int(run.data.metadata["PE_YIND"]) * int(meta["MYSUB"]) + np.arange(
            inner.sizes["y"]
        )
        selected = (global_y > meta["jyseps1_1"]) & (global_y <= meta["jyseps2_2"])
        core.append(inner.isel(y=np.flatnonzero(selected)))
        values.append((inner, outer, binner, bouter, selected))
    mean = xr.concat(core, dim="y").mean("y")
    result = []
    for run, (inner, outer, binner, bouter, selected) in zip(runs, values, strict=True):
        use_mean = run.data.options[MODEL].get_bool("phi_core_averagey")
        mask = xr.DataArray(selected, dims="y", coords={"y": inner.y})
        target = xr.where(mask, mean, inner) if use_mean else inner
        tau = (
            run.data.options[MODEL].evaluate_scalar("phi_boundary_timescale")
            * run.frequency
        )
        result.extend(
            [(binner, (target - binner) / tau), (bouter, (outer - bouter) / tau)]
        )
    return result


def ode_error(boundary, rate):
    dt = np.diff(boundary.t.values)
    require(np.allclose(dt, dt[0], rtol=1e-12, atol=0), "Irregular output cadence")
    axis = rate.get_axis_num("t")
    r = np.moveaxis(rate.values, axis, 0)
    cumulative = np.concatenate(
        [
            np.zeros_like(r[:1]),
            np.cumsum(
                (r[1:] + r[:-1]) * dt.reshape((-1,) + (1,) * (r.ndim - 1)) / 2, axis=0
            ),
        ]
    )
    b = np.moveaxis(boundary.values, boundary.get_axis_num("t"), 0)
    return float(np.max(abs(b - b[:1] - cumulative)))


def sheath_state(run):
    """Reconstruct the actual simple-sheath guards from saved interiors.

    This fixture explicitly uses exponential extrapolation (mode 1), no
    secondary emission, a zero wall potential and a single toroidal cell.
    """
    ds = run.data
    fields = {}
    for s in run.species:
        n = ds["N" + s.name]
        temperature = ds.Pe / ds.Ne if s.name == "e" else ds.Tion
        pressure = ds.Pe if s.name == "e" else n * temperature
        velocity = (
            sum(
                i.charge / i.mass * ds["NV" + i.name]
                for i in run.species
                if i.name != "e"
            )
            / ds.Ne
            if s.name == "e"
            else ds["NV" + s.name] / (s.mass * soft_floor(n, s.density_floor))
        )
        fields[s.name] = {
            k: v.copy(deep=True).transpose("t", "x", "y", "z")
            for k, v in {"n": n, "T": temperature, "p": pressure, "v": velocity}.items()
        }
    rank = int(ds.metadata["PE_YIND"])
    sources = {s.name: xr.zeros_like(ds.Eion) for s in run.species}
    total_heat = {s.name: xr.zeros_like(ds.Eion) for s in run.species}
    particle = {s.name: xr.zeros_like(ds.Eion) for s in run.species if s.name != "e"}
    physical_x = slice(run.mxg, -run.mxg)
    for sign, cell, neighbour, guard in (
        (-1, run.myg, run.myg + 1, run.myg - 1),
        (
            1,
            ds.sizes["y"] - run.myg - 1,
            ds.sizes["y"] - run.myg - 2,
            ds.sizes["y"] - run.myg,
        ),
    ):
        if (sign < 0 and rank != 0) or (sign > 0 and rank != 5):
            continue
        for s in run.species:
            f = fields[s.name]
            for name in ("n", "T", "p"):
                centre = f[name].isel(y=cell)
                next_cell = f[name].isel(y=neighbour)
                f[name].values[:, :, guard, :] = xr.where(
                    next_cell < 1e-10, centre, centre**2 / next_cell
                ).values
        electron = fields["e"]
        te = (electron["T"].isel(y=cell) + electron["T"].isel(y=guard)) / 2
        phi = 1.5 * ds.phi.isel(y=cell) - 0.5 * ds.phi.isel(y=neighbour)
        volume = ds.J * ds.dx * ds.dy * ds.dz
        for s in run.species:
            f = fields[s.name]
            n = (f["n"].isel(y=cell) + f["n"].isel(y=guard)) / 2
            temp = (f["T"].isel(y=cell) + f["T"].isel(y=guard)) / 2
            if s.name == "e":
                vface = (
                    sign
                    * np.sqrt(te / (2 * np.pi * s.mass))
                    * np.exp(-phi.clip(min=0) / te.clip(min=1e-5))
                )
                gamma = ds.options["sheath_boundary_simple"].evaluate_scalar("gamma_e")
            else:
                speed = np.sqrt(
                    (temp.clip(min=1e-5) + s.charge * te.clip(min=1e-5)) / s.mass
                )
                vface = (
                    np.maximum(f["v"].isel(y=cell), speed)
                    if sign > 0
                    else np.minimum(f["v"].isel(y=cell), -speed)
                )
                gamma = ds.options["sheath_boundary_simple"].evaluate_scalar("gamma_i")
            f["v"].values[:, :, guard, :] = (2 * vface - f["v"].isel(y=cell)).values
            # Preserve the current C++ geometry convention, including the ion
            # lower-target area evaluated with the neighbouring interior metric.
            other = neighbour if s.name != "e" and sign < 0 else guard
            area = (
                (ds.J.isel(y=cell) + ds.J.isel(y=other))
                / (np.sqrt(ds.g_22.isel(y=cell)) + np.sqrt(ds.g_22.isel(y=other)))
                * (ds.dx.isel(y=cell) + ds.dx.isel(y=other))
                / 2
                * (ds.dz.isel(y=cell) + ds.dz.isel(y=other))
                / 2
            )
            adv = (2.5 * temp + 0.5 * s.mass * vface**2) * n * vface
            heat = gamma * temp * n * vface
            for target, value in (
                (sources[s.name], -sign * (heat - adv) * area / volume.isel(y=cell)),
                (total_heat[s.name], -sign * heat * area / volume.isel(y=cell)),
            ):
                target.values[:, physical_x, cell, :] = value.transpose(
                    "t", "x", "z"
                ).values[:, physical_x, :]
            if s.name != "e":
                # Current diagnostics use minus the signed flux at both targets.
                value = -n * vface * area / volume.isel(y=cell)
                particle[s.name].values[:, physical_x, cell, :] = value.transpose(
                    "t", "x", "z"
                ).values[:, physical_x, :]
    return fields, sources, total_heat, particle


def audit(runs, manifest, diagnose):
    require(len(runs) == 6, "COMPASS fixture requires six ranks")
    tol = manifest["tolerances"]
    assertions = {}
    diagnostics = {}
    inventory = sum(energy_inventory(run) for run in runs)
    power_scale = runs[0].energy_scale * runs[0].length_scale ** 3 * runs[0].frequency
    power = lambda run, f: run.integrate(f) * power_scale
    sheath_power = xr.zeros_like(inventory.ion_saved)
    parallel_power = xr.zeros_like(sheath_power)
    parallel_faces = []
    raw_signal = 0.0
    for index, run in enumerate(runs):
        ds = run.data
        require(
            ds.metadata["NXPE"] == 1
            and ds.metadata["NYPE"] == 6
            and ds.sizes["z"] == 1,
            "Unsupported fixture decomposition",
        )
        require(
            ds.options[MODEL].get_bool("phi_boundary_relax")
            and not ds.options[MODEL].get_bool("poloidal_flows"),
            "Expected axisymmetric relaxation fixture",
        )
        densities = energy_densities(run)
        maximum = lambda f, run=run: float(abs(run.interior(f)).max())
        assertions[f"partition_{index}"] = (
            maximum(densities.ion_reconstructed - ds.Eion) < tol["field_absolute"]
        )
        assertions[f"forward_{index}"] = (
            maximum(densities.vorticity_reconstructed - ds.Vort) < tol["field_absolute"]
        )
        diagnostics[str(index)] = compare_diagnostics(run, densities=densities)
        for name, comparison in diagnostics[str(index)].items():
            assertions[f"diagnostic_{index}_{name}"] = comparison[
                "available"
            ] == diagnose and (
                not diagnose
                or comparison["max_absolute_error_normalised"] < tol["field_absolute"]
            )
        fields, sources, heat, particle = sheath_state(run)
        ions = [s for s in run.species if s.name != "e"]
        ion_source = sum(sources[s.name] for s in ions)
        sheath_power += power(run, sum(heat.values())).assign_coords(t=inventory.t)
        face = xr.zeros_like(ds.Eion)
        for s in ions:
            f = fields[s.name]
            energy = 2.5 * f["p"] + 0.5 * s.mass * f["n"] * f["v"] ** 2
            transport = parallel_flux(
                energy,
                f["v"],
                np.sqrt(f["T"] / s.mass),
                ds,
                myg=run.myg,
                periodic_y=True,
            )
            fp = transport.flux_y.copy(deep=True)
            area = ds.J * ds.dx * ds.dz / np.sqrt(ds.g_22)
            for cell in ([run.myg - 1] if index == 0 else []) + (
                [ds.sizes["y"] - run.myg - 1] if index == 5 else []
            ):
                value = (
                    (energy.isel(y=cell) + energy.isel(y=cell + 1))
                    / 2
                    * (f["v"].isel(y=cell) + f["v"].isel(y=cell + 1))
                    / 2
                    * (area.isel(y=cell) + area.isel(y=cell + 1))
                    / 2
                )
                fp.values[:, :, cell, :] = value.transpose("t", "x", "z").values
            face += fp
        parallel_faces.append(face)
        outward = (
            face.isel(y=ds.sizes["y"] - run.myg - 1)
            if index == 5
            else xr.zeros_like(face.isel(y=run.myg))
        ) - (
            face.isel(y=run.myg - 1)
            if index == 0
            else xr.zeros_like(face.isel(y=run.myg))
        )
        parallel_power += (
            outward.isel(x=slice(run.mxg, -run.mxg)).sum(("x", "z")) * power_scale
        ).assign_coords(t=inventory.t)
        if not diagnose:
            continue
        for s in run.species:
            assertions[f"sheath_{index}_{s.name}"] = (
                maximum(heat[s.name] - ds["E" + s.name + "_sheath"])
                < tol["field_absolute"]
            )
            if s.name != "e":
                assertions[f"particles_{index}_{s.name}"] = (
                    maximum(particle[s.name] - ds["S" + s.name + "_sheath"])
                    < tol["field_absolute"]
                )
        assertions[f"source_assembly_{index}"] = (
            maximum(ion_source - ds.Eion_src) < tol["field_absolute"]
        )
        chi = ds["ddt(Phi)"]
        coefficient = run.scalar_option("average_atomic_mass") / ds.Bxy**2
        raw = perpendicular_operator(coefficient, chi, ds, mxg=run.mxg)
        raw_signal = max(
            raw_signal,
            float(
                abs(
                    raw.flux_x.isel(
                        x=[run.mxg - 1, ds.sizes["x"] - run.mxg - 1],
                        y=slice(run.myg, -run.myg),
                    )
                ).max()
            ),
        )
        rho = sum(s.mass * ds["N" + s.name] for s in ions)
        weights = {
            "electric": ds.phi,
            "enthalpy": sum(
                s.mass / s.charge * 2.5 * ds["N" + s.name] * ds.Tion for s in ions
            )
            / rho,
            "kinetic": sum(
                s.mass
                / s.charge
                * 0.5
                * s.mass
                * ds["N" + s.name]
                * fields[s.name]["v"] ** 2
                for s in ions
            )
            / rho,
        }
        for s in ions:
            weights["particles_" + s.name] = s.mass / s.charge * ds["N" + s.name] / rho
            weights["momentum_" + s.name] = s.mass / s.charge * ds["NV" + s.name] / rho
        reconstructed = {}
        for name, weight in weights.items():
            op = perpendicular_operator(
                coefficient * weight, chi, ds, mxg=run.mxg, boundary_flux=False
            )
            assertions[f"closed_{index}_{name}"] = (
                float(
                    abs(
                        op.flux_x.isel(
                            x=[run.mxg - 1, ds.sizes["x"] - run.mxg - 1],
                            y=slice(run.myg, -run.myg),
                        )
                    ).max()
                )
                == 0
            )
            reconstructed[name] = -op.divergence
        assertions[f"closed_energy_{index}"] = (
            maximum(
                sum(reconstructed[n] for n in ("electric", "enthalpy", "kinetic"))
                - ds.DivEpol
            )
            < tol["field_absolute"]
        )
        closed = perpendicular_operator(
            coefficient, chi, ds, mxg=run.mxg, boundary_flux=False
        )
        assertions[f"closed_current_{index}"] = (
            maximum(-closed.divergence - ds.DivJpol) < tol["field_absolute"]
        )
        assertions[f"gauss_{index}"] = (
            float(abs(run.integrate(closed.divergence)).max())
            / float(run.integrate(xr.ones_like(ds.Eion)).max())
            < tol["field_absolute"]
        )
        assertions[f"axisymmetric_exb_{index}"] = (
            maximum(ds.DivEExB) < tol["field_absolute"]
            and maximum(ds.DivEReynolds) < tol["field_absolute"]
        )
    # Match actual topology connections instead of summing unrelated guards.
    sep = int(runs[0].data.metadata["ixseps1"])
    interfaces = [(i, i + 1, slice(sep, -runs[0].mxg)) for i in range(5)]
    interfaces += [(i, i + 1, slice(runs[0].mxg, sep)) for i in range(1, 4)]
    interfaces += [(4, 1, slice(runs[0].mxg, sep)), (0, 5, slice(runs[0].mxg, sep))]
    for a, b, x in interfaces:
        fa = parallel_faces[a].isel(x=x, y=runs[a].data.sizes["y"] - runs[a].myg - 1)
        fb = parallel_faces[b].isel(x=x, y=runs[b].myg - 1)
        assertions[f"interface_{a}_{b}_{x.start}"] = (
            float(abs(fa - fb).max()) < tol["field_absolute"]
        )
    ode = max(ode_error(b, r) for b, r in boundary_rates(runs))
    signal = max(float(abs(r).max()) for _, r in boundary_rates(runs))
    assertions["boundary_ode"] = ode < tol["boundary_ode_absolute"]
    assertions["relaxation_resolved"] = (
        signal > manifest["minimum_signals"]["boundary_relaxation"]
    )
    assertions["sheath_resolved"] = (
        float(abs(sheath_power).max()) > manifest["minimum_signals"]["sheath_power_W"]
    )
    report = {
        "assertions": assertions,
        "inventory_J": {n: inventory[n].values.tolist() for n in inventory},
        "diagnostics": diagnostics,
        "boundary_ode_max_normalised": ode,
        "relaxation_max_normalised": signal,
        "sheath_power_W": (-sheath_power).values.tolist(),
        "parallel_ion_outward_power_W": parallel_power.values.tolist(),
        "coverage": {
            "complete": False,
            "independent": [
                "inventory",
                "X-Z forward operator",
                "radial relaxation",
                "simple-sheath source assembly",
                "parallel ion boundary power and topology interfaces",
            ],
            "diagnostic_cross_checks": [
                "closed polarisation fluxes driven by saved auxiliary chi",
                "integrated ion RHS",
            ],
            "missing": [
                "general tokamak electrical/diamagnetic boundary work",
                "independent auxiliary-current RHS in shifted topology",
                "complete electron power/kinetic exchange audit",
            ],
        },
    }
    if diagnose:
        assertions["raw_current_resolved"] = (
            raw_signal > manifest["minimum_signals"]["raw_current"]
        )
        divparallel = sum(power(r, r.data.DivEparallel) for r in runs).assign_coords(
            t=inventory.t
        )
        reference = max(float(abs(divparallel).max()), 1e-30)
        report["parallel_power_relative_error"] = (
            float(abs(divparallel - parallel_power).max()) / reference
        )
        assertions["parallel_boundary_gauss"] = (
            report["parallel_power_relative_error"] < tol["power_normalised"]
        )
        rhs = sum(
            power(
                r,
                r.data.Eion_src
                - r.data.phi_Div_Je
                - sum(
                    r.data[n]
                    for n in (
                        "DivEExB",
                        "DivEReynolds",
                        "DivEelectric",
                        "DivEparallel",
                        "DivEpol",
                    )
                ),
            )
            for r in runs
        ).assign_coords(t=inventory.t)
        rhs.attrs["units"] = "W"
        budget = energy_balance(
            inventory,
            reservoir="ion_saved",
            source_power=rhs,
            outward_power=xr.zeros_like(rhs).assign_attrs(units="W"),
        )
        budget.attrs["complete"] = False
        budget.attrs["coverage"] = (
            "diagnostic RHS identity; incomplete independent tokamak boundary/source audit"
        )
        report["diagnostic_ion_balance"] = histories(budget)
        assertions["diagnostic_rhs_balance"] = (
            float(abs(budget.relative_residual).max()) < tol["balance_relative"]
        )
    report["passed"] = all(assertions.values())
    return report, inventory


def main(inputs):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=ROOT / "build/hermes-3")
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--profile", choices=("smoke", "ci"), default="ci")
    args = parser.parse_args()
    manifest = json.loads((inputs / "criteria.json").read_text())
    grid = inputs / "data/compass-36x48.grd.nc"
    if not grid.is_file():
        grid = ROOT / "examples/tokamak-2D/compass-36x48.grd.nc"
    require(
        grid.is_file() and args.exe.is_file(),
        "Missing COMPASS grid or built executable",
    )
    fixture = inputs / "data/BOUT.inp"
    options = BoutOptionsFile(str(fixture), gridfilename=str(grid))
    require(
        {n.strip() for n in options["hermes"]["components"].split(",")}
        == {"e", "d+", "impurity", MODEL, "sheath_boundary_simple"},
        "Unaccounted components",
    )
    for name in ("Pe", "Nd+", "Nimpurity"):
        require(
            options[name].evaluate_scalar("source") == 0,
            "Fixture requires explicit zero imposed sources",
        )
    output = Path(tempfile.mkdtemp(prefix=inputs.name + "-", dir=args.output_root))
    print(f"Results: {output}", flush=True)
    summary = {"passed": False, "runs": [], "manifest": manifest}
    write_json(output / "result.json", summary)
    cadences = (
        manifest["cadences_normalised"][-2:]
        if args.profile == "smoke"
        else manifest["cadences_normalised"]
    )
    enabled = None
    previous = None
    previous_inventory = None
    refinement = []
    for cadence, diagnose in [*((c, True) for c in cadences), (cadences[-1], False)]:
        folder = output / f"dt{cadence}-diagnose-{'on' if diagnose else 'off'}"
        print(f"Running {folder.name} on six ranks", flush=True)
        execute(
            args.exe.resolve(),
            folder,
            fixture,
            manifest["case"],
            cadence,
            manifest["duration_normalised"],
            diagnose,
            args.timeout,
            ranks=manifest["ranks"],
            files=(grid,),
        )
        runs = read_ranks(folder)
        report, inventory = audit(runs, manifest, diagnose)
        write_json(folder / "energy-report.json", report)
        summary["runs"].append(str(folder))
        write_json(output / "result.json", summary)
        require(
            report["passed"],
            f"Failed {[n for n, v in report['assertions'].items() if not v]}: {folder / 'energy-report.json'}",
        )
        if diagnose:
            residual = max(
                abs(v) for v in report["diagnostic_ion_balance"]["relative_residual"]
            )
            if previous is not None:
                ratio = residual / previous
                refinement.append(ratio)
                require(
                    ratio < manifest["tolerances"]["maximum_refinement_ratio"],
                    "Sampled ion RHS balance did not converge with cadence",
                )
                xr.testing.assert_allclose(
                    previous_inventory,
                    inventory.sel(t=previous_inventory.t),
                    rtol=1e-7,
                    atol=1e-10,
                )
            previous = residual
            previous_inventory = inventory
        if cadence == cadences[-1] and diagnose:
            enabled = runs
        if not diagnose:
            for on, off in zip(enabled, runs, strict=True):
                for name in ESSENTIAL:
                    xr.testing.assert_allclose(
                        on.data[name],
                        off.data[name],
                        rtol=manifest["tolerances"]["history_relative"],
                        atol=manifest["tolerances"]["history_absolute"],
                    )
            stripped = []
            for on in enabled:
                data = on.data.drop_vars(
                    [
                        n
                        for n in on.data
                        if "t" in on.data[n].dims
                        and n not in ESSENTIAL
                        and n != "t_array"
                    ]
                )
                stripped.append(
                    EnergyRun(
                        data,
                        on.species,
                        on.model,
                        on.mxg,
                        on.myg,
                        on.energy_scale,
                        on.length_scale,
                        on.frequency,
                        on.source,
                    )
                )
            stripped_report, stripped_inventory = audit(stripped, manifest, False)
            require(stripped_report["passed"], "Stripped-diagnostic audit failed")
            xr.testing.assert_allclose(stripped_inventory, inventory)
    summary.update(
        passed=True,
        diagnostics_independent=True,
        optional_fields_stripped=True,
        independent_tokamak_power_coverage_complete=False,
        rhs_balance_refinement_ratios=refinement,
    )
    write_json(output / "result.json", summary)
    print(
        f"{inputs.name} passed (explicitly incomplete general tokamak power coverage)"
    )
    return 0


def run_test(inputs):
    try:
        return main(inputs)
    except (ValueError, AssertionError) as error:
        print(f"Tokamak energy test failed: {error}", file=sys.stderr)
        return 1
