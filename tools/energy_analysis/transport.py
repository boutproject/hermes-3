"""Discrete slab transport, independent of Hermes and BOUT operator calls.

Fluxes are normalised, integrated over a whole cell face, oriented towards
increasing coordinate, and indexed by the cell on the lower side of the face.
Outer guard results are undefined. Complete halos are required. Polarisation
weights use zero-gradient physical wall guards. The orchestration layer checks supported metrics.
"""

import numpy as np
import xarray as xr

from .dataset import AuditError
from .operators import (
    face_divergence,
    face_work,  # noqa: F401 -- public transport API
    neighbor,
    perpendicular_faces,
    weighted_faces,  # noqa: F401 -- public transport API
)


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
    field,
    stream,
    geometry,
    *,
    mxg=2,
    boundary_flux=False,
    weight=1,
    limiter="mc",
    periodic_x=False,
):
    """Hermes X-Z MC ExB transport, with weights applied after limiting.

    By default matches the shared firstX/lastX switch even for periodic X.
    Reynolds transport requests periodic_x=True to retain seam flow.
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
    if periodic_x:
        pass  # Complete periodic faces, including the firstX/lastX seam.
    elif boundary_flux:
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


def reynolds_face_fluxes(
    c,
    vorticity,
    phi,
    psi,
    geometry,
    *,
    mxg=2,
    periodic_x=False,
    limiter="mc",
):
    """Full-vorticity MC advection plus the central metric remainder R.

    Potential and metric guards must already be valid. Tangential diagonal
    divergences use the original potential samples, including wall jumps.
    No derived guards or artificial wall-current closure are imposed.
    """
    xp, xm = lambda f: neighbor(f, "x"), lambda f: neighbor(f, "x", -1)
    zp, zm = lambda f: neighbor(f, "z"), lambda f: neighbor(f, "z", -1)
    avg = lambda f: (f + xp(f) + zp(f) + xp(zp(f))) / 4
    derivative_x = lambda f: (xp(f) - xm(f)) / (2 * geometry.dx)
    derivative_z = lambda f: (zp(f) - zm(f)) / (2 * geometry.dz)
    Phi = phi + psi
    c = xr.ones_like(Phi) * c
    A = c * geometry.J
    g11, g13, g33 = geometry.g11, geometry.g13, geometry.g33
    D = g11 * g33 - g13**2
    Ax, Az, Dx, Dz = derivative_x(A), derivative_z(A), derivative_x(D), derivative_z(D)
    a, b = derivative_x(g11), derivative_z(g11)
    d, e = derivative_x(g13), derivative_z(g13)
    f, g = derivative_x(g33), derivative_z(g33)
    coefficients = (
        -(A * b + Az * g11) / 2,
        -A * a / 2 - A * e - Ax * g11 / 2 - Az * g13,
        (
            A * D * a * g33
            + 2 * A * D * g13 * g
            - A * Dx * g11 * g33
            - A * Dz * g13 * g33
            + Az * D * g13 * g33
        )
        / (2 * D * g33),
        (A * D * g - A * Dx * g13 - A * Dz * g33 - Ax * D * g13) / (2 * D),
        -(A * D * a - A * Dx * g11 - A * Dz * g13 - Az * D * g13) / (2 * D),
        -(
            A * D * g11 * g
            + 2 * A * D * a * g13
            - A * Dx * g11 * g13
            - A * Dz * g11 * g33
            + Ax * D * g11 * g13
        )
        / (2 * D * g11),
        A * d + A * g / 2 + Ax * g13 + Az * g33 / 2,
        (A * f + Ax * g33) / 2,
    )
    corner_phi, corner_Phi = avg(phi), avg(Phi)
    hx_corner = (xp(Phi) - Phi + xp(zp(Phi)) - zp(Phi)) / (2 * avg(geometry.dx))
    hz_corner = (zp(Phi) - Phi + xp(zp(Phi)) - xp(Phi)) / (2 * avg(geometry.dz))
    # Diagonal fluxes exclude cross-metric terms by construction.
    gx = xr.Dataset({**{k: geometry[k] for k in geometry}, "g13": xr.zeros_like(g13)})
    diagonal = perpendicular_faces(c, Phi, gx, mxg=mxg)
    volume = geometry.J * geometry.dx * geometry.dy * geometry.dz
    Lx = (diagonal["x"] - xm(diagonal["x"])) / volume
    Lz = (diagonal["z"] - zm(diagonal["z"])) / volume
    fluxes = {}
    for dim, shift, tangential, offset in (("x", xp, Lz, 0), ("z", zp, Lx, 4)):
        face = lambda f, shift=shift: (f + shift(f)) / 2
        dx, dz = face(geometry.dx) if dim == "x" else geometry.dx, face(geometry.dz)
        J = face(geometry.J) if dim == "x" else geometry.J
        if dim == "x":
            velocity = J * (corner_phi - zm(corner_phi))
            px, pz = (xp(phi) - phi) / dx, velocity / (J * dz)
            hx, hz = (xp(Phi) - Phi) / dx, (corner_Phi - zm(corner_Phi)) / dz
            mixed = (hx_corner - zm(hx_corner)) / dz
            remainder = -J * face(c) * (face(g11) * px + face(g13) * pz) * mixed
            remainder -= J * face(g13) / face(g33) * px * face(tangential)
            area = dz
        else:
            velocity = J * (xm(corner_phi) - corner_phi)
            px, pz = -velocity / (J * dx), (zp(phi) - phi) / dz
            hx, hz = (corner_Phi - xm(corner_Phi)) / dx, (zp(Phi) - Phi) / dz
            mixed = (hz_corner - xm(hz_corner)) / dx
            remainder = J * face(c) * (face(g13) * px + face(g33) * pz) * mixed
            remainder += J * face(g13) / face(g11) * pz * face(tangential)
            area = dx
        remainder += sum(
            face(coefficients[offset + i]) * p
            for i, p in enumerate((hx * px, hx * pz, hz * px, hz * pz))
        )
        if limiter == "mc":
            left, right = mc_edges(vorticity, dim)
            donor = xr.where(velocity > 0, right, shift(left))
        elif limiter == "central":
            # Operator-only convergence probe; production retains MC donors.
            donor = face(vorticity)
        else:
            raise ValueError(f"Unknown Reynolds limiter: {limiter}")
        fluxes[dim] = (
            velocity * (donor - face(tangential)) + area * remainder
        ) * geometry.dy
    return fluxes


def fv_xz(
    coefficient, potential, geometry, *, mxg=2, boundary_flux=True, periodic_x=False
):
    return flux_dataset(
        perpendicular_faces(
            coefficient,
            potential,
            geometry,
            mxg=mxg,
            boundary_flux=boundary_flux,
            periodic_x=periodic_x,
        ),
        geometry,
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
