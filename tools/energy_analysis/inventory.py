"""Independent ion/electron reservoirs and optional diagnostic comparisons."""

import numpy as np
import xarray as xr

from .dataset import AuditError
from .operators import perpendicular_operator


def soft_floor(value, floor):
    positive = value.clip(min=0)
    return positive + floor * np.exp(-positive / floor)


def _kinetic(run, species):
    ds = run.data
    n = ds["N" + species.name]
    momentum = "NV" + species.name
    if momentum in ds:
        if ds[momentum].attrs.get("units_type") != "hermes":
            raise AuditError(f"Unknown normalisation for {momentum}")
        velocity = ds[momentum] / (species.mass * soft_floor(n, species.density_floor))
        return 0.5 * species.mass * n * velocity**2
    if "V" + species.name in ds:
        # Fixed/prescribed flows can expose velocity instead of momentum.
        velocity = ds["V" + species.name]
        if velocity.attrs.get("units_type") != "hermes":
            raise AuditError(f"Unknown normalisation for V{species.name}")
        return 0.5 * species.mass * n * velocity**2
    if any(
        t in species.components
        for t in (
            "evolve_momentum",
            "zero_current",
            "electromagnetic",
            "fixed_velocity",
        )
    ):
        raise AuditError(
            f"Missing saved momentum/velocity for flowing species {species.name}"
        )
    return xr.zeros_like(n)


def energy_densities(run):
    """Return normalised reservoirs; no optional diagnostics are accessed."""
    ds = run.data
    ions = [s for s in run.species if s.name != "e" and abs(s.charge) > 1e-5]
    abar = run.scalar_option("average_atomic_mass")
    psi = sum(s.mass / s.charge * ds["N" + s.name] * ds.Tion for s in ions) / abar
    coefficient = abar / ds.Bxy**2
    operators = perpendicular_operator(coefficient, ds.phi + psi, ds, mxg=run.mxg)
    ion_parallel = sum(_kinetic(run, s) for s in ions)
    ion_thermal = 1.5 * sum(ds["N" + s.name] * ds.Tion for s in ions)
    electron = next(s for s in run.species if s.name == "e")
    result = xr.Dataset(
        {
            "ion_saved": ds.Eion,
            "ion_thermal": ion_thermal,
            "ion_parallel": ion_parallel,
            "ion_perpendicular": operators.kperp,
            "ion_reconstructed": ion_thermal + ion_parallel + operators.kperp,
            "electron_thermal": 1.5 * ds.Pe,
            "electron_parallel": _kinetic(run, electron)
            if "evolve_momentum" in electron.components
            else xr.zeros_like(ds.Pe),
            "pressure_potential": psi,
            "vorticity_reconstructed": operators.divergence,
        }
    )
    result["plasma_saved"] = (
        result.ion_saved + result.electron_thermal + result.electron_parallel
    )
    result["plasma_reconstructed"] = (
        result.ion_reconstructed + result.electron_thermal + result.electron_parallel
    )
    for name in result:
        result[name].attrs = {
            "units_type": "hermes",
            "provenance": "saved" if name == "ion_saved" else "reconstructed",
        }
    return result


def energy_inventory(run, *, densities=None):
    """Volume-integrated reservoirs in joules and time in seconds."""
    densities = energy_densities(run) if densities is None else densities
    names = [
        n
        for n in densities
        if n not in ("pressure_potential", "vorticity_reconstructed")
    ]
    scale = run.energy_scale * run.length_scale**3
    result = xr.Dataset(
        {name: run.integrate(densities[name]) * scale for name in names}
    )
    result = result.assign_coords(t=run.data.t / run.frequency)
    result.t.attrs = {"units": "s"}
    for name in result:
        result[name].attrs = {
            "units": "J",
            "provenance": densities[name].attrs["provenance"],
        }
    result.attrs = {
        "source": run.source,
        "coverage": "energy inventory; flux/source reconstruction not yet implemented",
    }
    return result


def compare_diagnostics(run, *, densities=None):
    """Physical-cell and volume-integrated differences for available diagnostics.

    Absence is recorded, never interpreted as a zero diagnostic. Errors are
    normalised field/volume errors, so iteration tolerances can be applied.
    """
    densities = energy_densities(run) if densities is None else densities
    comparisons = {
        "ion_KE_perp": "ion_perpendicular",
        "ion_KE_parallel": "ion_parallel",
        "ion_pressure_potential": "pressure_potential",
    }
    result = {}
    for diagnostic, reconstructed in comparisons.items():
        if diagnostic not in run.data:
            result[diagnostic] = {"available": False}
            continue
        error = densities[reconstructed] - run.data[diagnostic]
        absolute = float(abs(run.interior(error)).max().compute())
        scale = float(abs(run.interior(run.data[diagnostic])).max().compute())
        result[diagnostic] = {
            "available": True,
            "max_absolute_error_normalised": absolute,
            "max_relative_error": absolute / max(scale, np.finfo(float).tiny),
            "max_integral_error_normalised": float(
                abs(run.integrate(error)).max().compute()
            ),
        }
    return result
