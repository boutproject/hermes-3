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


def energy_balance(inventory, *, source_power=None, outward_power=None):
    """R = delta E - integral Q + integral P, in joules.

    Supplied powers must be scalar time series in W on exactly this time axis.
    Unknown contributions leave the residual unavailable. Explicit zero powers
    are appropriate only for a case with independently known zero contributions.
    """
    dt = regular_timestep(inventory.t)
    energy = inventory.plasma_saved
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
    for name in result:
        if name != "relative_residual":
            result[name].attrs["units"] = "J"
    result.attrs = {
        "complete": complete,
        "coverage": "supplied source and boundary powers"
        if complete
        else "source/boundary coverage incomplete",
        "time_spacing_s": dt,
    }
    return result
