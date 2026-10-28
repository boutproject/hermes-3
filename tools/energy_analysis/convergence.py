"""Coordinate-aware X-Z restriction and self-convergence diagnostics."""

import numpy as np
from scipy.interpolate import CubicSpline

from .dataset import AuditError


def overlap_matrix(source, target, *, source_width, target_width, period=None):
    """Conservative cell-average map, including shifted periodic cell faces.

    Rows preserve constants; width-weighted columns preserve the integral.
    Both grids must cover the same interval (or one complete periodic box).
    """
    source, target = np.asarray(source), np.asarray(target)
    sw, tw = float(source_width), float(target_width)
    if sw <= 0 or tw <= 0 or not np.isfinite(np.r_[source, target, sw, tw]).all():
        raise AuditError("Invalid cell coordinates/widths")
    if np.any(np.diff(source) <= 0) or np.any(np.diff(target) <= 0):
        raise AuditError("Cell coordinates must increase")
    if not np.isclose(len(source) * sw, len(target) * tw, rtol=1e-12):
        raise AuditError("Restriction grids cover different lengths")
    if period is None and not np.isclose(
        source[0] - sw / 2, target[0] - tw / 2, atol=tw * 1e-10
    ):
        raise AuditError("Nonperiodic restriction grids have different origins")
    if period is not None and not np.isclose(period, len(source) * sw, rtol=1e-12):
        raise AuditError("Periodic restriction must cover one whole period")
    result = np.zeros((len(target), len(source)))
    for shift in [-period, 0, period] if period is not None else [0]:
        low = np.maximum(target[:, None] - tw / 2, source[None, :] + shift - sw / 2)
        high = np.minimum(target[:, None] + tw / 2, source[None, :] + shift + sw / 2)
        result += np.maximum(high - low, 0) / tw
    if not np.allclose(result.sum(1), 1, atol=1e-11) or not np.allclose(
        result.sum(0) * tw, sw, atol=sw * 1e-10
    ):
        raise AuditError("Restriction failed constant/conservation checks")
    return result


def restrict(
    field,
    source_x,
    source_z,
    target_x,
    target_z,
    *,
    lx,
    lz,
    representation="cell_average",
):
    """Map a physical X-Z array without index-based spatial shifts."""
    f = np.asarray(field)
    if f.shape != (len(source_x), len(source_z)) or not np.isfinite(f).all():
        raise AuditError("Invalid reference field")
    if representation == "point":
        # Fourier interpolation in periodic Z, followed by cubic interpolation X.
        n = len(source_z)
        wave = np.fft.fftfreq(n, d=lz / n) * 2 * np.pi
        phases = np.exp(
            1j * wave[:, None] * (np.asarray(target_z) - source_z[0])[None, :]
        )
        sampled_z = (np.fft.fft(f, axis=1) @ phases / n).real
        return CubicSpline(source_x, sampled_z, axis=0)(target_x)
    if representation != "cell_average":
        raise AuditError(f"Unknown field representation {representation}")
    rx = overlap_matrix(
        source_x,
        target_x,
        source_width=lx / len(source_x),
        target_width=lx / len(target_x),
    )
    rz = overlap_matrix(
        source_z,
        target_z,
        source_width=lz / len(source_z),
        target_width=lz / len(target_z),
        period=lz,
    )
    return rx @ f @ rz.T


def error_norms(field, reference, volume):
    error, ref, weight = (
        np.asarray(field) - reference,
        np.asarray(reference),
        np.asarray(volume),
    )
    if (
        not np.isfinite(error).all()
        or not np.isfinite(weight).all()
        or np.any(weight <= 0)
    ):
        raise AuditError("Error norms require finite fields and positive volumes")
    l2 = float(np.sqrt(np.sum(error**2 * weight) / np.sum(weight)))
    linf = float(np.max(abs(error)))
    rms = float(np.sqrt(np.sum(ref**2 * weight) / np.sum(weight)))
    maximum = float(np.max(abs(ref)))
    return {
        "L2": l2,
        "Linf": linf,
        "relative_L2": l2 / rms if rms else None,
        "relative_Linf": linf / maximum if maximum else None,
    }


def global_differences(values, refinement):
    q = np.asarray(values, dtype=float)
    if refinement < 2 or not np.isfinite(q).all() or len(q) < 2:
        raise AuditError(
            "Global convergence requires >=2 finite grids and refinement>=2"
        )
    signed = q[:-1] - q[1:]
    floor = 128 * np.finfo(float).eps * max(np.max(abs(q)), np.finfo(float).tiny)
    rows = []
    for i, value in enumerate(q):
        order, reason = (
            None,
            "finest reference"
            if i == len(q) - 1
            else "two adjacent differences required",
        )
        if i < len(signed) - 1:
            if abs(signed[i]) <= floor or abs(signed[i + 1]) <= floor:
                reason = "roundoff dominated"
            elif signed[i] * signed[i + 1] <= 0:
                reason = "oscillatory differences"
            else:
                order = float(
                    np.log(abs(signed[i] / signed[i + 1])) / np.log(refinement)
                )
                reason = "adjacent signed differences"
        rows.append(
            {
                "Q_J": float(value),
                "adjacent_signed_J": float(signed[i]) if i < len(signed) else None,
                "finest_difference_J": float(abs(value - q[-1]))
                if i < len(signed)
                else None,
                "observed_order": order,
                "order_status": reason,
            }
        )
    return rows
