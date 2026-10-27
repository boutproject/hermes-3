"""Compact machine-readable inventory reports."""

import subprocess
import sys
from dataclasses import asdict
from importlib.metadata import version
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


def inventory_report(run):
    densities = energy_densities(run)
    inventory = energy_inventory(run, densities=densities).compute()
    dt = regular_timestep(inventory.t)
    partition = run.interior(densities.ion_reconstructed - run.data.Eion)
    forward = run.interior(densities.vorticity_reconstructed - run.data.Vort)
    import xhermes

    checkout = Path(xhermes.__file__).resolve().parents[1]
    report = {
        "schema_version": 1,
        "source": run.source,
        "python": sys.executable,
        "xhermes": {
            "version": version("xhermes"),
            "path": str(checkout),
            **revision(checkout),
        },
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
        "partition_max_absolute_error_normalised": float(
            abs(partition).max().compute()
        ),
        "partition_error_by_time_normalised": abs(partition)
        .max(["x", "y", "z"])
        .compute()
        .values.tolist(),
        "forward_max_absolute_error_normalised": float(abs(forward).max().compute()),
        "forward_error_by_time_normalised": abs(forward)
        .max(["x", "y", "z"])
        .compute()
        .values.tolist(),
        "diagnostics": compare_diagnostics(run, densities=densities),
        "balance": dict(energy_balance(inventory).attrs),
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
    families = ("DivEExB", "DivEReynolds", "DivEelectric", "DivEparallel", "DivEpol")
    boundary = sum(flux_power(run, fluxes, n) for n in families).assign_attrs(
        units="W", provenance="reconstructed ion face powers"
    )
    source = (
        integrated_power(
            run,
            (fluxes.Eion_src - fluxes.phi_Div_Je).assign_attrs(
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
