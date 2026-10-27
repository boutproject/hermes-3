"""Regular-output cumulative energy balances with explicit term coverage."""

import numpy as np
import xarray as xr

from .dataset import AuditError


def regular_timestep(time, *, rtol=1e-9):
    """Validate finite, increasing, regularly spaced saved times.

    Spacing tolerance combines relative interval error and 32 roundoff units at
    the absolute time scale, to allow subtraction of large restart timestamps.
    """
    values = np.asarray(time, dtype=float)
    if values.ndim != 1 or values.size < 2 or not np.isfinite(values).all():
        raise AuditError("Time axis must contain at least two finite samples")
    intervals = np.diff(values)
    if np.any(intervals <= 0):
        raise AuditError(
            "Time axis must be strictly increasing; check restart duplicates"
        )
    dt = float(np.median(intervals))
    atol = 32 * np.finfo(float).eps * max(abs(values).max(), dt)
    if not np.allclose(intervals, dt, rtol=rtol, atol=atol):
        raise AuditError("Time axis is irregular or has output gaps")
    return dt


def energy_balance(
    inventory, *, source_power=None, outward_power=None, reservoir="plasma_saved"
):
    """R = delta E - integral Q + integral P, in joules.

    Select a saved/reconstructed reservoir by name; plasma_saved is the default.
    Supplied powers must be scalar time series in W on exactly this time axis.
    Unknown contributions leave the residual unavailable. Explicit zero powers
    are appropriate only for a case with independently known zero contributions.
    """
    dt = regular_timestep(inventory.t)
    energy = inventory[reservoir]
    if (
        energy.dims != ("t",)
        or energy.attrs.get("units") != "J"
        or not np.isfinite(energy.values).all()
    ):
        raise AuditError("Energy must be a finite scalar time series in J")
    change = energy - energy.isel(t=0)
    result = xr.Dataset({"energy_change": change})
    complete = source_power is not None and outward_power is not None
    for name, power in (
        ("source_energy", source_power),
        ("outward_energy", outward_power),
    ):
        if power is None:
            continue
        if power.dims != ("t",) or not power.t.equals(inventory.t):
            raise AuditError(
                "Power must have the inventory time axis and only dimension t"
            )
        if power.attrs.get("units") != "W" or not np.isfinite(power.values).all():
            raise AuditError("Power must be finite and explicitly expressed in W")
        values = power.values
        result[name] = ("t", np.r_[0, np.cumsum((values[1:] + values[:-1]) * dt / 2)])
    if complete:
        result["residual"] = change - result.source_energy + result.outward_energy
        scale = max(float(abs(energy).max()), np.finfo(float).tiny)
        result["relative_residual"] = result.residual / scale
        result.relative_residual.attrs["units"] = "1"
        result["interval_residual"] = result.residual.diff("t").reindex(t=inventory.t)
        result["interval_relative_residual"] = result.interval_residual / scale
        result.interval_relative_residual.attrs["units"] = "1"
    for name in result:
        if name not in ("relative_residual", "interval_relative_residual"):
            result[name].attrs["units"] = "J"
    result.attrs = {
        "complete": complete,
        "coverage": "supplied source and boundary powers"
        if complete
        else "source/boundary coverage incomplete",
        "time_spacing_s": dt,
        "reservoir": reservoir,
    }
    if complete:
        result.attrs["reference_energy_J"] = scale
        result.attrs["source_provenance"] = source_power.attrs.get(
            "provenance", "supplied"
        )
        result.attrs["boundary_provenance"] = outward_power.attrs.get(
            "provenance", "supplied"
        )
    return result


def cadence_comparison(coarse, fine):
    """Compare cumulative residuals on nested output axes at fixed duration.

    A second-order trapezoid estimate of the fine history's sampling error is
    |R_coarse-R_fine|/(ratio**2-1). A persistent common residual is thereby
    visible separately from the cadence-dependent part. This estimate assumes
    smooth power histories; callers must verify convergence with another run.
    """
    coarse_dt, fine_dt = (regular_timestep(ds.t) for ds in (coarse, fine))
    ratio = coarse_dt / fine_dt
    stride = round(ratio)
    if stride < 2 or not np.isclose(ratio, stride, rtol=1e-9):
        raise AuditError("Cadences must have an integer refinement ratio >= 2")
    sampled = fine.isel(t=slice(None, None, stride))
    time_atol = fine_dt * 1e-9 + 32 * np.finfo(float).eps * max(
        float(abs(coarse.t).max()), float(abs(fine.t).max())
    )
    if sampled.sizes["t"] != coarse.sizes["t"] or not np.allclose(
        sampled.t, coarse.t, rtol=0, atol=time_atol
    ):
        raise AuditError("Cadence runs must share start, end and nested saved times")
    for ds in (coarse, fine):
        if not ds.attrs.get("complete") or "residual" not in ds:
            raise AuditError("Cadence comparison requires complete balances")
        if (
            ds.residual.attrs.get("units") != "J"
            or not np.isfinite(ds.residual.values).all()
        ):
            raise AuditError("Cadence residuals must be finite and expressed in J")
    if coarse.attrs.get("reservoir") != fine.attrs.get("reservoir"):
        raise AuditError("Cadence comparison requires the same energy reservoir")
    difference = np.asarray(coarse.residual) - np.asarray(sampled.residual)
    extrapolated = np.asarray(sampled.residual) - difference / (stride**2 - 1)
    return {
        "refinement_ratio": stride,
        "max_residual_difference_J": float(abs(difference).max()),
        "estimated_fine_sampling_error_J": float(
            abs(difference).max() / (stride**2 - 1)
        ),
        "coarse_max_residual_J": float(abs(coarse.residual).max()),
        "fine_max_residual_J": float(abs(fine.residual).max()),
        "fine_final_residual_J": float(fine.residual.isel(t=-1)),
        "estimated_cadence_independent_max_residual_J": float(abs(extrapolated).max()),
        "provenance": "nested-output trapezoid refinement; smooth-power assumption",
    }
