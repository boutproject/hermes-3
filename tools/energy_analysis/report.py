"""Compact machine-readable inventory reports."""

import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .balance import energy_balance, regular_timestep
from .inventory import compare_diagnostics, energy_densities, energy_inventory


def revision(path):
    def git(*args):
        result = subprocess.run(
            ["git", "-C", str(path), *args], capture_output=True, text=True, check=False
        )
        return result.stdout.strip() if result.returncode == 0 else None

    return {
        "revision": git("rev-parse", "HEAD"),
        "dirty": bool(git("status", "--porcelain")),
    }


def summary(run):
    """JSON-serialisable inventories, recovery, drift and power coverage."""
    from scipy.integrate import cumulative_trapezoid

    from .budgets import closed_balance, particle_inventory, source_balance
    from .dataset import AuditError
    from .sources import external_power
    from .transfers import (
        collisional_heat_exchange,
        electron_curvature_work,
        integrated_power,
    )
    from .walls import wall_budget

    densities = energy_densities(run)
    inventory = energy_inventory(run, densities=densities).compute()
    particles = particle_inventory(run).compute()
    drift, power, walls = closed_balance(run), external_power(run), wall_budget(run)
    report = {
        "source": run.source,
        "time_s": inventory.t.values.tolist(),
        "inventory_J": {name: inventory[name].values.tolist() for name in inventory},
        "particles": {name: particles[name].values.tolist() for name in particles},
        "closed_balance": {
            "coverage": drift.attrs["coverage"],
            "missing_terms": drift.attrs["missing_terms"],
            "maxima": {
                name: float(abs(drift[name]).max())
                for name in drift
                if name.endswith("relative_drift")
            },
        },
        "external_power": {
            "coverage": power.attrs["coverage"],
            "missing_terms": power.attrs["missing_terms"],
            "terms": {name: power[name].values.tolist() for name in power},
        },
        "wall_budget": {
            "coverage": walls.attrs["coverage"],
            "missing_terms": walls.attrs["missing_terms"],
        },
        "diagnostics": compare_diagnostics(run, densities=densities),
    }
    if run.formulation == "total_energy":
        report["recovery_error"] = float(
            abs(run.interior(densities.ion_reconstructed - run.data.Eion)).max()
        )
    for name in ("ion_perpendicular", "ion_parallel"):
        report[name + "_change_relative"] = float(
            abs(inventory[name] - inventory[name].isel(t=0)).max()
            / max(abs(float(inventory.plasma_native.isel(t=0))), np.finfo(float).tiny)
        )
    try:
        exchange = integrated_power(run, electron_curvature_work(run))
        report["curvature_exchange_relative"] = float(
            np.max(abs(cumulative_trapezoid(exchange, x=inventory.t, initial=0)))
        ) / max(abs(float(inventory.plasma_native.isel(t=0))), np.finfo(float).tiny)
    except AuditError as error:
        report["curvature_exchange"] = {
            "coverage": "incomplete",
            "missing_terms": [str(error)],
        }
    exchange = collisional_heat_exchange(run)
    if exchange:
        report["collisional_heat_exchange"] = {
            "coverage": exchange.attrs["coverage"],
            "missing_terms": exchange.attrs["missing_terms"],
            "power_W": {name: exchange[name].values.tolist() for name in exchange},
        }
    if inventory.sizes["t"] >= 2:
        balance = source_balance(run)
        report["source_balance"] = {
            "coverage": balance.attrs["coverage"],
            "missing_terms": balance.attrs["missing_terms"],
            "maxima": {
                name: float(abs(balance[name]).max())
                for name in balance
                if name.endswith("relative_residual")
            },
        }
    return report


def inventory_report(run):
    densities = energy_densities(run)
    inventory = energy_inventory(run, densities=densities).compute()
    dt = regular_timestep(inventory.t)
    partition = (
        run.interior(densities.ion_reconstructed - run.data.Eion)
        if run.formulation == "total_energy"
        else None
    )
    forward = run.interior(densities.vorticity_reconstructed - run.data.Vort)
    report = {
        "schema_version": 1,
        "source": run.source,
        "python": sys.executable,
        "hermes_revision": run.data.metadata.get("HERMES_REVISION"),
        "analysis_checkout": revision(Path(__file__).resolve().parents[2]),
        "run_id": run.data.metadata.get("run_id"),
        "normalisations": {
            "energy_density_J_m3": run.energy_scale,
            "length_m": run.length_scale,
            "frequency_s_inv": run.frequency,
        },
        "domain": {"mxg": run.mxg, "myg": run.myg, "dimensions": dict(run.data.sizes)},
        "species": [asdict(s) for s in run.species],
        "time_s": inventory.t.values.tolist(),
        "time_spacing_s": dt,
        "inventory_J": {name: inventory[name].values.tolist() for name in inventory},
        "partition_max_absolute_error_normalised": float(abs(partition).max().compute())
        if partition is not None
        else None,
        "partition_error_by_time_normalised": abs(partition)
        .max(["x", "y", "z"])
        .compute()
        .values.tolist()
        if partition is not None
        else None,
        "forward_max_absolute_error_normalised": float(abs(forward).max().compute()),
        "forward_error_by_time_normalised": abs(forward)
        .max(["x", "y", "z"])
        .compute()
        .values.tolist(),
        "diagnostics": compare_diagnostics(run, densities=densities),
        "balance": dict(energy_balance(inventory, reservoir="plasma_native").attrs),
    }
    if not all(np.isfinite(inventory[name]).all() for name in inventory):
        from .dataset import AuditError

        raise AuditError("Nonfinite reconstructed inventory; check saved guard fields")
    return report, inventory


def flux_report(run, *, sources=None):
    """Compact independent flux/balance report; unsupported coverage stays explicit.

    Sources are an explicit PrescribedSources manifest, never inferred from
    absent diagnostics. Detailed face maps remain available through
    reconstruct_fluxes; this report retains only scalar histories and errors.
    """
    from .dataset import AuditError
    from .fluxes import (
        compare_flux_diagnostics,
        flux_power,
        reconstruct_fluxes,
        reconstruct_sources,
    )
    from .transfers import integrated_power

    def histories(ds):
        return {
            **ds.attrs,
            **{
                name: [float(v) if np.isfinite(v) else None for v in ds[name].values]
                for name in ds
            },
        }

    report, inventory = inventory_report(run)
    try:
        fluxes = reconstruct_fluxes(run, sources=sources)
    except AuditError as error:
        report["flux_audit"] = {
            "complete": False,
            "reason": str(error),
            "provenance": "unavailable",
        }
        return report
    families = (
        "DivEExB",
        "DivEelectric",
        "DivEparallel",
        "DivEpol",
        "DivEviscous",
    )
    boundary = sum(flux_power(run, fluxes, n) for n in families).assign_attrs(
        units="W", provenance="reconstructed ion face powers"
    )
    source = (
        integrated_power(
            run,
            (
                fluxes.Eion_src - fluxes.phi_Div_Je - fluxes.electron_viscous_heating
            ).assign_attrs(
                units_type="hermes",
                provenance="prescribed external sources plus reconstructed electron exchange",
            ),
        )
        if sources
        else None
    )
    ion = energy_balance(
        inventory, reservoir="ion_saved", source_power=source, outward_power=boundary
    )
    details = {
        "coverage": fluxes.attrs,
        "ion_balance": histories(ion),
        "diagnostics": compare_flux_diagnostics(run, fluxes),
        "outward_ion_power_W": boundary.values.tolist(),
        "reynolds_power_W": {
            name: integrated_power(run, fluxes[name]).values.tolist()
            for name in ("reynolds_work", "reynolds_heating")
        },
        "gauss_error_W": {
            n: float(
                abs(integrated_power(run, fluxes[n]) - flux_power(run, fluxes, n))
                .max()
                .compute()
            )
            for n in families
        },
    }
    if sources:
        external = reconstruct_sources(run, sources)
        electron = next(s for s in run.species if s.name == "e")
        outward = (
            boundary
            + flux_power(run, fluxes, "electron_ExB")
            + flux_power(run, fluxes, "electron_parallel_thermal")
        )
        reason = None
        if "evolve_momentum" in electron.components:
            if not fluxes.attrs["periodic_y"]:
                reason = "Open electron kinetic boundary flux is unavailable"
            elif (
                float(
                    abs(run.interior(run.data.phi - run.data.phi.isel(x=run.mxg, z=0)))
                    .max()
                    .compute()
                )
                > 1e-12
            ):
                reason = "Electron kinetic ExB flux reconstruction is unavailable"
            else:
                outward += flux_power(run, fluxes, "electron_parallel_kinetic")
        plasma = energy_balance(
            inventory,
            source_power=integrated_power(run, external.plasma_src),
            outward_power=outward.assign_attrs(
                units="W", provenance="reconstructed ion and electron faces"
            )
            if reason is None
            else None,
        )
        details["plasma_balance"] = histories(plasma)
        if reason:
            details["plasma_balance"]["reason"] = reason
        if reason is None:
            instantaneous = (
                integrated_power(
                    run,
                    (
                        fluxes["ddt(Eion)"]
                        + 1.5 * fluxes["ddt(Pe)"]
                        + fluxes.electron_kinetic_rhs
                        - external.plasma_src
                    ).assign_attrs(units_type="hermes"),
                )
                + outward
            )
            details["plasma_instantaneous_unexplained_power_W"] = (
                instantaneous.values.tolist()
            )
    details["complete"] = ion.attrs["complete"]
    report["flux_audit"] = details
    return report
