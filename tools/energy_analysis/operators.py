"""X-Z finite-volume operators in Hermes normalised coordinates.

Each flux is integrated over the transverse coordinate in one Y plane, points
towards increasing X or Z, and is indexed by the cell on its lower side.
Outer guard results are undefined (NaN); callers must select physical cells.
Z is periodic. X requires at least one saved guard on each side.
"""

import numpy as np
import xarray as xr

from .dataset import AuditError


def spectral_ddz(field, dz):
    """Periodic FFT first derivative using the saved coordinate spacing.

    BOUT's FFT first derivative has zero response at the even-grid Nyquist
    mode. Z must be uniform; complete Z lines are rechunked while other
    dimensions remain lazy. This is distinct from the finite-volume Kperp.
    """
    if "z" not in field.dims:
        raise AuditError("FFT derivative requires a named z dimension")
    spacing = np.asarray(dz)
    if (
        not np.isfinite(spacing).all()
        or np.any(spacing <= 0)
        or not np.allclose(spacing, spacing.flat[0], rtol=1e-12, atol=0)
    ):
        raise AuditError("FFT derivative requires uniform positive saved dz")
    n = field.sizes["z"]
    wave = 2 * np.pi * np.fft.rfftfreq(n, d=float(spacing.flat[0]))
    if n % 2 == 0:
        wave[-1] = 0

    def derivative(values):
        return np.fft.irfft(1j * wave * np.fft.rfft(values, axis=-1), n=n, axis=-1)

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


def neighbor(field, dim, direction=1):
    if dim not in field.dims or (dim == "y" and field.sizes[dim] == 1):
        return field
    if dim == "z":
        return field.roll(z=-direction, roll_coords=False)
    return field.shift({dim: -direction})


def face_divergence(fluxes, geometry):
    """Discrete Gauss divergence of area-integrated coordinate face fluxes."""
    volume = geometry.J * geometry.dx * geometry.get("dy", 1) * geometry.dz
    result = sum(f - neighbor(f, dim, -1) for dim, f in fluxes.items()) / volume
    result.attrs = {"units_type": "hermes", "provenance": "reconstructed"}
    return result


def face_work(fluxes, potential, geometry):
    """Half the sum of oriented face fluxes times potential jumps per volume."""
    work = 0
    for dim, flux in fluxes.items():
        jump = neighbor(potential, dim) - potential
        work += flux * jump + neighbor(flux * jump, dim, -1)
    return work / (2 * geometry.J * geometry.dx * geometry.get("dy", 1) * geometry.dz)


def close_radial_faces(faces, field, mxg):
    """Close the two physical wall rows after evaluating the geometric faces."""
    ix = xr.DataArray(np.arange(field.sizes["x"]), dims="x", coords={"x": field.x})
    return {
        **faces,
        "x": faces["x"].where((ix != mxg - 1) & (ix != field.sizes["x"] - mxg - 1), 0),
    }


def neumann_weights(weight, *, mxg=2, periodic_x=False):
    """Local copy with interior values on every physical radial guard column."""
    if periodic_x:
        return weight
    w = weight.copy()
    for cell, step in ((mxg, -1), (w.sizes["x"] - mxg - 1, 1)):
        for layer in range(1, mxg + 1):
            w[{"x": cell + step * layer}] = w.isel(x=cell)
    return w


def weighted_faces(faces, weight, *, mxg=2, periodic_x=False):
    """Face-average a Neumann wall weight on the common polarisation current."""
    w = neumann_weights(weight, mxg=mxg, periodic_x=periodic_x)
    return {d: f * (w + neighbor(w, d)) / 2 for d, f in faces.items()}


def perpendicular_faces(
    coefficient, potential, geometry, *, mxg=2, boundary_flux=True, periodic_x=False
):
    """Area-integrated diagonal/cross fluxes, computed once per face and corner."""
    if not {"x", "z"}.issubset(potential.dims):
        raise ValueError("X-Z operator requires named x and z dimensions")
    if mxg < 1 or potential.sizes["x"] <= 2 * mxg:
        raise ValueError("X-Z operator requires retained radial guards")
    a = xr.ones_like(potential) * coefficient
    f = potential
    J, dx, dz = (geometry[k] for k in ("J", "dx", "dz"))
    g11, g33, g13 = (geometry[k] for k in ("g11", "g33", "g13"))
    xp = lambda v: v.shift(x=-1) if "x" in v.dims else v
    xm = lambda v: v.shift(x=1) if "x" in v.dims else v
    zp = lambda v: v.roll(z=-1, roll_coords=False) if "z" in v.dims else v
    zm = lambda v: v.roll(z=1, roll_coords=False) if "z" in v.dims else v
    corner = (
        (a + xp(a) + zp(a) + xp(zp(a)))
        * (J * g13 + xp(J * g13) + zp(J * g13) + xp(zp(J * g13)))
        / 16
    )
    cross_x = corner * (zp(f) - f + xp(zp(f)) - xp(f)) / 2
    cross_z = corner * (xp(f) - f + xp(zp(f)) - zp(f)) / 2
    fx = (a + xp(a)) / 2 * (J * g11 * dz + xp(J * g11 * dz)) * (xp(f) - f) / (
        dx + xp(dx)
    ) + (cross_x + zm(cross_x)) / 2
    fz = (a + zp(a)) / 2 * (J * g33 * dx + zp(J * g33 * dx)) * (zp(f) - f) / (
        dz + zp(dz)
    ) + (cross_z + xm(cross_z)) / 2
    if not boundary_flux and not periodic_x:
        index = xr.DataArray(range(f.sizes["x"]), dims="x", coords={"x": f.x})
        fx = fx.where((index != mxg - 1) & (index != f.sizes["x"] - mxg - 1), 0)
    return {"x": fx * geometry.get("dy", 1), "z": fz * geometry.get("dy", 1)}


def perpendicular_operator(
    coefficient, potential, geometry, *, mxg=2, boundary_flux=True, periodic_x=False
):
    """Divergence, geometric Kperp and shared area-integrated X-Z faces."""
    faces = perpendicular_faces(coefficient, potential, geometry, mxg=mxg)
    ke = face_work(faces, potential, geometry) / 2
    if not boundary_flux and not periodic_x:
        faces = close_radial_faces(faces, potential, mxg)
    return xr.Dataset(
        {
            "divergence": face_divergence(faces, geometry),
            "kperp": ke,
            **{"flux_" + d: f / geometry.get("dy", 1) for d, f in faces.items()},
        },
        attrs={"units_type": "hermes", "flux_orientation": "increasing coordinate"},
    )


def calculate_kperp(coefficient, potential, geometry, *, mxg=2):
    """Scheme-matched perpendicular energy density (normalised)."""
    return (
        face_work(
            perpendicular_faces(coefficient, potential, geometry, mxg=mxg),
            potential,
            geometry,
        )
        / 2
    )
