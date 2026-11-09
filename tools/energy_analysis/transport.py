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
    if dim not in field.dims or (dim == "y" and field.sizes[dim] == 1):
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


def xppm_faces(
    field, stream, geometry, *, mxg=2, boundary_flux=False, weight=1, limiter="mc"
):
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
    if limiter == "mc":
        left_x, right_x = mc_edges(field, "x")
        left_z, right_z = mc_edges(field, "z")
    elif limiter is None:
        # Kinetic advection uses central faces; tests also apply them to the current.
        left_x = right_x = (field + xp(field)) / 2
        left_z = right_z = (field + zp(field)) / 2
    else:
        raise ValueError(f"Unknown XPPM limiter: {limiter}")
    fx = vx * (right_x if limiter is None else xr.where(vx > 0, right_x, xp(left_x)))
    fz = vz * (right_z if limiter is None else xr.where(vz > 0, right_z, zp(left_z)))
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
    return {"x": fx, "z": fz}


def xppm(field, stream, geometry, *, mxg=2, boundary_flux=False, weight=1):
    """Divergence and faces of MC XPPM; weights follow reconstruction."""
    return flux_dataset(
        xppm_faces(
            field, stream, geometry, mxg=mxg, boundary_flux=boundary_flux, weight=weight
        ),
        geometry,
    )


def central_advection_faces(field, stream, geometry, *, mxg=2):
    """Central face averages times the XPPM corner-stream velocities."""
    return xppm_faces(field, stream, geometry, mxg=mxg, limiter=None)


def face_work(fluxes, potential, geometry):
    """Half the sum of oriented face fluxes times potential jumps per volume."""
    work = 0
    for dim, flux in fluxes.items():
        jump = neighbor(potential, dim) - potential
        work += flux * jump + neighbor(flux * jump, dim, -1)
    return work / (2 * geometry.J * geometry.dx * geometry.dy * geometry.dz)


def solenoidal_reynolds_faces(c, phi, Phi, geometry):
    """Corner curl of -c J grad(phi).grad(Phi)/2, including cross metric."""
    xp, zp = lambda f: neighbor(f, "x"), lambda f: neighbor(f, "z")
    avg = lambda f: (f + xp(f) + zp(f) + xp(zp(f))) / 4
    gradient = lambda f: (
        (xp(f) - f + xp(zp(f)) - zp(f)) / (2 * avg(geometry.dx)),
        (zp(f) - f + xp(zp(f)) - xp(f)) / (2 * avg(geometry.dz)),
    )
    px, pz = gradient(phi)
    hx, hz = gradient(Phi)
    stream = (
        -avg(c * geometry.J)
        / 2
        * (
            avg(geometry.g11) * px * hx
            + avg(geometry.g13) * (px * hz + pz * hx)
            + avg(geometry.g33) * pz * hz
        )
    )
    return {
        "x": (stream - neighbor(stream, "z", -1)) * geometry.dy,
        "z": -(stream - neighbor(stream, "x", -1)) * geometry.dy,
    }


def reynolds_face_fluxes(
    c,
    vorticity,
    phi,
    psi,
    geometry,
    *,
    mxg=2,
    walls=True,
    periodic_x=False,
    limiter="mc",
    fill_guards=None,
):
    """Skew X-Z Reynolds faces; closed split wall flux retains the corner curl.

    The orchestration layer supplies communication/Y guard filling; standalone
    tests can supply analytic periodic halos. Physical X scalar guards follow
    the component's midpoint-zero Q and linear-extrapolated Lambda rules.
    """
    Phi = phi + psi
    D = geometry.g11 * geometry.g33 - geometry.g13**2
    lam = (
        c
        * D
        * geometry.J
        / 2
        * perpendicular_operator(
            1 / (geometry.J * D),
            phi,
            geometry,
            mxg=mxg,
            boundary_flux=True,
            periodic_x=periodic_x,
        ).divergence
    )
    Q = arakawa(phi, psi, geometry.dx, geometry.dz)
    if fill_guards is not None:
        Q, lam = fill_guards(Q), fill_guards(lam)
    if not periodic_x:
        # Copy before writes: broadcast arrays and dask-backed values are allowed.
        Q, lam = Q.copy(), lam.copy()
        for cell, step in ((mxg, -1), (phi.sizes["x"] - mxg - 1, 1)):
            for layer in range(1, mxg + 1):
                Q[{"x": cell + step * layer}] = -Q.isel(x=cell)
                lam[{"x": cell + step * layer}] = lam.isel(x=cell) + layer * (
                    lam.isel(x=cell) - lam.isel(x=cell - step)
                )
    omega = xppm_faces(vorticity / 2, phi, geometry, mxg=mxg, limiter=limiter)
    rotation = xppm_faces(lam, Phi, geometry, mxg=mxg, limiter=limiter)
    q = fv_xz(c / 2, Q, geometry, mxg=mxg, boundary_flux=False, periodic_x=periodic_x)
    curl = solenoidal_reynolds_faces(c, phi, Phi, geometry)
    split = {d: omega[d] + rotation[d] + q["flux_" + d] for d in ("x", "z")}
    if walls and not periodic_x:
        ix = xr.DataArray(np.arange(phi.sizes["x"]), dims="x", coords={"x": phi.x})
        split["x"] = split["x"].where(
            (ix != mxg - 1) & (ix != phi.sizes["x"] - mxg - 1), 0
        )
    return {d: split[d] + curl[d] for d in ("x", "z")}


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
        0.5
        * (ar + vr)
        * nr
        * (velocity - yp(velocity))
        * (vr - (velocity + yp(velocity)) / 2)
        * (area + yp(area))
        / 2
    )
    left = (
        0.5
        * (al - vl)
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
