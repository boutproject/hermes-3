"""Diagnostic-independent conservation and source-driven balance entry points."""

import numpy as np
import xarray as xr
from scipy.integrate import cumulative_trapezoid

from .balance import energy_balance
from .inventory import energy_inventory
from .sources import coverage, external_power
from .walls import wall_budget


def particle_inventory(run):
    values = {
        "particles_" + s.name: run.integrate(run.data["N" + s.name])
        * (run.data.metadata["Nnorm"] * run.length_scale**3)
        for s in run.species
        if "N" + s.name in run.data
    }
    values["charge"] = sum(
        s.charge * values["particles_" + s.name]
        for s in run.species
        if "particles_" + s.name in values
    )
    values["vorticity"] = run.integrate(run.data.Vort)
    return coverage(xr.Dataset(values).assign_coords(t=run.data.t / run.frequency))


def closed_balance(run):
    """Relative drifts from record zero; coverage verifies a closed system."""
    inventory = energy_inventory(run)
    particles = particle_inventory(run)
    values = {"energy": inventory["plasma_native"], **dict(particles.data_vars)}
    result = xr.Dataset()
    for name, value in values.items():
        drift = value - value.isel(t=0)
        scale = max(abs(float(value.isel(t=0))), np.finfo(float).tiny)
        if name in ("charge", "vorticity"):
            scale = (
                max(
                    float(
                        sum(
                            abs(particles["particles_" + s.name].isel(t=0))
                            for s in run.species
                            if "particles_" + s.name in particles
                        )
                    ),
                    np.finfo(float).tiny,
                )
                if name == "charge"
                else max(float(abs(run.integrate(abs(run.data.Vort))).max()), 1.0)
            )
        result[name + "_change"] = drift
        result[name + "_relative_drift"] = drift / scale
    power, walls = external_power(run), wall_budget(run)
    missing = power.attrs["missing_terms"] + walls.attrs["missing_terms"]
    if float(abs(power.power).max()) > 0 or any(
        float(abs(power[n]).max()) > 0 for n in power if n.startswith("particles_")
    ):
        missing.append("nonzero external sources")
    if float(abs(walls.outward).max()) > 0:
        missing.append("open wall power")
    return coverage(result, missing)


def source_balance(run, reservoir="plasma"):
    """Cumulative delta inventory minus integrated external powers and wall losses.

    An incomplete balance retains inventories and known power histories, but
    omits residuals; unknown contributions are never treated as zero.
    """
    inventory, power, walls = (
        energy_inventory(run),
        external_power(run),
        wall_budget(run),
    )
    missing = power.attrs["missing_terms"] + walls.attrs["missing_terms"]
    if reservoir == "plasma":
        reservoir = "plasma_native"
    result = energy_balance(
        inventory,
        reservoir=reservoir,
        source_power=power.power if not missing else None,
        outward_power=walls.outward if not missing else None,
    )
    result["known_external_power"] = power.power
    result["known_outward_power"] = walls.outward
    if not missing:
        particles = particle_inventory(run)
        for s in run.species:
            name = "particles_" + s.name
            if name not in particles or name not in power:
                continue
            rates = power[name] - walls.get(name, xr.zeros_like(power[name]))
            supplied = xr.DataArray(
                cumulative_trapezoid(rates, x=particles.t, initial=0),
                dims="t",
                coords={"t": particles.t},
            )
            result[name + "_residual"] = (
                particles[name] - particles[name].isel(t=0) - supplied
            )
            result[name + "_supplied"] = supplied
            result[name + "_relative_residual"] = result[name + "_residual"] / max(
                float(abs(particles[name]).max()), np.finfo(float).tiny
            )
    return coverage(result, missing)
