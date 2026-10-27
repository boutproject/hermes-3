"""X-Z finite-volume operators in Hermes normalised coordinates.

Each flux is integrated over the transverse coordinate in one Y plane, points
towards increasing X or Z, and is indexed by the cell on its lower side.
Outer guard results are undefined (NaN); callers must select physical cells.
Z is periodic. X requires at least one saved guard on each side.
"""

import xarray as xr


def perpendicular_operator(
    coefficient, potential, geometry, *, mxg=2, boundary_flux=True, periodic_x=False
):
    """Return divergence, kinetic energy and shared diagonal/cross face fluxes.

    Mixed terms use one symmetric bilinear form at each X-Z corner, including
    spatially varying coefficients and metrics. Calculations remain lazy.
    """
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
    cell_area = J * dx * dz
    # Kperp always uses the unmodified geometric flux, as in the C++ operator.
    ke = (
        (xp(f) - f) * fx
        + (f - xm(f)) * xm(fx)
        + (zp(f) - f) * fz
        + (f - zm(f)) * zm(fz)
    ) / (4 * cell_area)
    if not boundary_flux and not periodic_x:
        index = xr.DataArray(range(f.sizes["x"]), dims="x", coords={"x": f.x})
        fx = fx.where((index != mxg - 1) & (index != f.sizes["x"] - mxg - 1), 0)
    div = (fx - xm(fx) + fz - zm(fz)) / cell_area
    return xr.Dataset(
        {"divergence": div, "kperp": ke, "flux_x": fx, "flux_z": fz},
        attrs={"units_type": "hermes", "flux_orientation": "increasing coordinate"},
    )


def calculate_kperp(coefficient, potential, geometry, *, mxg=2):
    """Scheme-matched perpendicular energy density (normalised)."""
    return perpendicular_operator(coefficient, potential, geometry, mxg=mxg)["kperp"]
