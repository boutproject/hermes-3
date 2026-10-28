"""Discrete slab transport, independent of Hermes and BOUT operator calls.

Fluxes are normalised, integrated over a whole cell face, oriented towards
increasing coordinate, and indexed by the cell on the lower side of the face.
Outer guard results are undefined. Complete halos are required; no physical
boundary values are invented. The orchestration layer checks supported metrics.
"""

import numpy as np
import xarray as xr

from .dataset import AuditError
from .operators import perpendicular_operator


def neighbor(field, dim, direction=1):
    if dim not in field.dims:
        return field
    if dim == "z":
        return field.roll(z=-direction, roll_coords=False)
    return field.shift({dim: -direction})


def mc_edges(field, dim):
    """Monotonized-central reconstruction, including zero slope at extrema."""
    plus, minus = neighbor(field, dim), neighbor(field, dim, -1)
    candidates = (2 * (plus - field), (plus - minus) / 2, 2 * (field - minus))
    positive = (candidates[0] > 0) & (candidates[1] > 0) & (candidates[2] > 0)
    negative = (candidates[0] < 0) & (candidates[1] < 0) & (candidates[2] < 0)
    magnitude = np.minimum(
        np.minimum(abs(candidates[0]), abs(candidates[1])), abs(candidates[2])
    )
    slope = xr.where(positive, magnitude, xr.where(negative, -magnitude, 0))
    return field - slope / 2, field + slope / 2


def face_divergence(fluxes, geometry):
    """Discrete Gauss divergence of area-integrated coordinate face fluxes."""
    volume = geometry.J * geometry.dx * geometry.dy * geometry.dz
    result = sum(f - neighbor(f, dim, -1) for dim, f in fluxes.items()) / volume
    result.attrs = {"units_type": "hermes", "provenance": "reconstructed"}
    return result


def flux_dataset(fluxes, geometry):
    return xr.Dataset(
        {
            "divergence": face_divergence(fluxes, geometry),
            **{"flux_" + d: f for d, f in fluxes.items()},
        },
        attrs={
            "units_type": "hermes",
            "flux_orientation": "increasing coordinate; lower-cell index",
            "flux_measure": "area integrated",
        },
    )


def xppm(field, stream, geometry, *, mxg=2, boundary_flux=False, weight=1):
    """Hermes X-Z MC ExB transport, with weights applied after limiting.

    Matches the current firstX/lastX boundary switch even for periodic X.
    Uniform geometry in X-Z is required by the orchestration layer, because
    the C++ Z scatter uses the donor cell's volume for both neighbors.
    """
    if mxg < 2:
        raise AuditError("XPPM requires two retained radial halos")
    xp = lambda v: neighbor(v, "x")
    zp = lambda v: neighbor(v, "z")
    zm = lambda v: neighbor(v, "z", -1)
    corner = (stream + xp(stream) + zp(stream) + xp(zp(stream))) / 4
    vx = (geometry.J + xp(geometry.J)) / 2 * (corner - zm(corner)) / geometry.dz
    vz = geometry.J * (neighbor(corner, "x", -1) - corner) / geometry.dx
    left_x, right_x = mc_edges(field, "x")
    left_z, right_z = mc_edges(field, "z")
    fx = vx * xr.where(vx > 0, right_x, xp(left_x))
    fz = vz * xr.where(vz > 0, right_z, zp(left_z))
    ix = xr.DataArray(np.arange(field.sizes["x"]), dims="x", coords={"x": field.x})
    low, high = mxg - 1, field.sizes["x"] - mxg - 1
    if boundary_flux:
        midpoint = (field + xp(field)) / 2
        fx = xr.where(ix == low, vx * xr.where(vx < 0, xp(left_x), midpoint), fx)
        fx = xr.where(ix == high, vx * xr.where(vx > 0, right_x, midpoint), fx)
    else:
        fx = fx.where((ix != low) & (ix != high), 0)
    w = xr.ones_like(field) * weight
    fx *= (w + xp(w)) / 2 * geometry.dy * geometry.dz
    fz *= (w + zp(w)) / 2 * geometry.dx * geometry.dy
    return flux_dataset({"x": fx, "z": fz}, geometry)


def fv_xz(
    coefficient, potential, geometry, *, mxg=2, boundary_flux=True, periodic_x=False
):
    result = perpendicular_operator(
        coefficient,
        potential,
        geometry,
        mxg=mxg,
        boundary_flux=boundary_flux,
        periodic_x=periodic_x,
    )
    return flux_dataset(
        {d: result["flux_" + d] * geometry.dy for d in ("x", "z")}, geometry
    )


def arakawa(f, g, dx, dz):
    xp, xm = lambda v: neighbor(v, "x"), lambda v: neighbor(v, "x", -1)
    zp, zm = lambda v: neighbor(v, "z"), lambda v: neighbor(v, "z", -1)
    jpp = (zp(f) - zm(f)) * (xp(g) - xm(g)) - (xp(f) - xm(f)) * (zp(g) - zm(g))
    jpx = (
        xp(g) * (xp(zp(f)) - xp(zm(f)))
        - xm(g) * (xm(zp(f)) - xm(zm(f)))
        - zp(g) * (xp(zp(f)) - xm(zp(f)))
        + zm(g) * (xp(zm(f)) - xm(zm(f)))
    )
    jxp = (
        xp(zp(g)) * (zp(f) - xp(f))
        - xm(zm(g)) * (xm(f) - zm(f))
        - xm(zp(g)) * (zp(f) - xm(f))
        + xp(zm(g)) * (xp(f) - zm(f))
    )
    return (jpp + jpx + jxp) / (12 * dx * dz)


def spectral_d2z(field, dz):
    spacing = np.asarray(dz)
    if (
        not np.isfinite(spacing).all()
        or np.any(spacing <= 0)
        or not np.allclose(spacing, spacing.flat[0], rtol=1e-12, atol=0)
    ):
        raise AuditError("Spectral second derivative requires uniform positive dz")
    n = field.sizes["z"]
    wave = 2 * np.pi * np.fft.rfftfreq(n, d=float(spacing.flat[0]))

    def derivative(values):
        return np.fft.irfft(-(wave**2) * np.fft.rfft(values, axis=-1), n=n, axis=-1)

    if field.chunks is not None:
        field = field.chunk({"z": -1})
    return xr.apply_ufunc(
        derivative,
        field,
        input_core_dims=[["z"]],
        output_core_dims=[["z"]],
        dask="parallelized",
        output_dtypes=[float],
    ).transpose(*field.dims)


def central_parallel(field, geometry):
    """Uniform-B identity-coordinate Div_par/Grad_par, centred C2 DDY."""
    if field.sizes["y"] == 1:
        return xr.zeros_like(field)
    return (neighbor(field, "y") - neighbor(field, "y", -1)) / (
        2 * geometry.dy * np.sqrt(geometry.g_22)
    )


def parallel_flux(
    field,
    velocity,
    wave,
    geometry,
    *,
    myg=2,
    periodic_y=True,
    fixflux=True,
    momentum=False,
):
    """MC split parallel transport; optionally Div_par_fvv (density argument).

    Physical energy boundaries use midpoint f and velocity when fixflux=True.
    The momentum boundary variant is intentionally unsupported here.
    """
    if myg < 2:
        if field.sizes["y"] == 1:
            return flux_dataset({"y": xr.zeros_like(field)}, geometry)
        raise AuditError("Parallel MC fluxes require two retained Y halos")
    if momentum and not periodic_y:
        raise AuditError("Open momentum boundary flux reconstruction is unavailable")
    yp = lambda v: neighbor(v, "y")
    fl, fr = mc_edges(field, "y")
    vl, vr = mc_edges(velocity, "y")
    speed = np.maximum(
        np.maximum(wave, yp(wave)), np.maximum(abs(velocity), abs(yp(velocity)))
    )
    flux = 0.5 * fr * (vr + speed) * (vr if momentum else 1) + 0.5 * yp(fl) * (
        yp(vl) - speed
    ) * (yp(vl) if momentum else 1)
    if not periodic_y:
        iy = xr.DataArray(np.arange(field.sizes["y"]), dims="y", coords={"y": field.y})
        low, high = myg - 1, field.sizes["y"] - myg - 1
        midpoint = (field + yp(field)) / 2
        vface = (velocity + yp(velocity)) / 2
        if fixflux:
            boundary = midpoint * vface
            flux = xr.where((iy == low) | (iy == high), boundary, flux)
        else:
            flux = xr.where(
                iy == low, yp(fl) * vface - yp(wave) * (yp(fl) - midpoint), flux
            )
            flux = xr.where(iy == high, fr * vface + wave * (fr - midpoint), flux)
    area = geometry.J * geometry.dx * geometry.dz / np.sqrt(geometry.g_22)
    return flux_dataset({"y": flux * (area + yp(area)) / 2}, geometry)


def numerical_parallel_heating(density, velocity, wave, geometry):
    """Div_par_fvv_heating interior expression; multiply by species mass."""
    ym, yp = lambda v: neighbor(v, "y", -1), lambda v: neighbor(v, "y")
    nl, nr = mc_edges(density, "y")
    vl, vr = mc_edges(velocity, "y")
    ar = np.maximum(
        np.maximum(wave, yp(wave)), np.maximum(abs(velocity), abs(yp(velocity)))
    )
    al = np.maximum(
        np.maximum(wave, ym(wave)), np.maximum(abs(velocity), abs(ym(velocity)))
    )
    area = geometry.J * geometry.dx * geometry.dz / np.sqrt(geometry.g_22)
    right = (
        (ar + 0.5 * vr)
        * nr
        * (velocity - yp(velocity))
        * (vr - (velocity + yp(velocity)) / 2)
        * (area + yp(area))
        / 2
    )
    left = (
        (al - 0.5 * vl)
        * nl
        * (velocity - ym(velocity))
        * (vl - (velocity + ym(velocity)) / 2)
        * (area + ym(area))
        / 2
    )
    return (right + left) / (geometry.J * geometry.dx * geometry.dy * geometry.dz)


def parallel_kinetic_flux(density, velocity, wave, geometry, *, myg=2):
    """Kinetic flux matching MC momentum, density and numerical heating.

    Interior/periodic Y only. Each split contribution uses the receiving
    cell's velocity, as in BOUT's viscous-heating flux diagnostic.
    """
    if myg < 2:
        raise AuditError("Kinetic parallel flux requires two Y halos")
    yp = lambda f: neighbor(f, "y")
    nl, nr = mc_edges(density, "y")
    vl, vr = mc_edges(velocity, "y")
    a = np.maximum(
        np.maximum(wave, yp(wave)), np.maximum(abs(velocity), abs(yp(velocity)))
    )
    flux = 0.5 * nr * (vr + a) * (
        yp(velocity) * vr - 0.5 * yp(velocity) ** 2
    ) + 0.5 * yp(nl) * (yp(vl) - a) * (velocity * yp(vl) - 0.5 * velocity**2)
    area = geometry.J * geometry.dx * geometry.dz / np.sqrt(geometry.g_22)
    return flux_dataset({"y": flux * (area + yp(area)) / 2}, geometry)


def numerical_heating_defect(density, velocity, wave, geometry):
    """Excess heating in BOUT's current Div_par_fvv_heating implementation.

    The split particle flux is n_edge*(a +/- v_edge)/2, whereas heating uses
    n_edge*(a +/- v_edge/2). Their difference is n_edge*a/2. This is a tracked
    discretisation contribution, not a prescribed external source.
    """
    ym, yp = lambda f: neighbor(f, "y", -1), lambda f: neighbor(f, "y")
    nl, nr = mc_edges(density, "y")
    vl, vr = mc_edges(velocity, "y")
    ar = np.maximum(
        np.maximum(wave, yp(wave)), np.maximum(abs(velocity), abs(yp(velocity)))
    )
    al = np.maximum(
        np.maximum(wave, ym(wave)), np.maximum(abs(velocity), abs(ym(velocity)))
    )
    area = geometry.J * geometry.dx * geometry.dz / np.sqrt(geometry.g_22)
    right = (
        0.5
        * ar
        * nr
        * (velocity - yp(velocity))
        * (vr - (velocity + yp(velocity)) / 2)
        * (area + yp(area))
        / 2
    )
    left = (
        0.5
        * al
        * nl
        * (velocity - ym(velocity))
        * (vl - (velocity + ym(velocity)) / 2)
        * (area + ym(area))
        / 2
    )
    return (right + left) / (geometry.J * geometry.dx * geometry.dy * geometry.dz)
