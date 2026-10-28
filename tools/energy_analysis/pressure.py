"""Native FD-X/FFT-Z energy and legacy viscous work in an orthogonal slab.

The Z contribution is the local quadratic form -Phi*D2Z(Phi)/2. Its integral
includes the Nyquist mode, unlike the square of BOUT's FFT first derivative.
Radial half-face contributions use saved physical guards. This density is a
representation of the native quadratic form, not a positive local reservoir.
"""

import numpy as np
import xarray as xr

from .dataset import AuditError
from .transport import neighbor, spectral_d2z


def legacy_operator(coefficient, potential, geometry, *, mxg=2):
    for name in ("J", "dx", "dz", "g11", "g33", "Bxy"):
        values = np.asarray(geometry[name])
        if not np.isfinite(values).all() or not np.allclose(
            values, values.flat[0], rtol=1e-12
        ):
            raise AuditError(f"Legacy native energy requires uniform {name}")
    if float(abs(geometry.g13).max()) > 1e-14:
        raise AuditError("Legacy native energy currently requires g13=0")
    a = np.asarray(coefficient)
    if not np.allclose(a, a.flat[0], rtol=1e-12):
        raise AuditError("Legacy native energy requires a constant coefficient")
    f, dx = potential, geometry.dx
    xp, xm = neighbor(f, "x"), neighbor(f, "x", -1)
    d2z = spectral_d2z(f, geometry.dz)
    kx = coefficient * geometry.g11 * ((xp - f) ** 2 + (f - xm) ** 2) / (4 * dx**2)
    kz = -0.5 * coefficient * geometry.g33 * f * d2z
    return xr.Dataset(
        {
            "kperp": kx + kz,
            "divergence": coefficient
            * (geometry.g11 * (xp - 2 * f + xm) / dx**2 + geometry.g33 * d2z),
        }
    )


def legacy_viscous_work(run):
    """Actual corner-gradient heating and its native kinetic-energy work.

    Their sum exposes the heating truncation defect. Boundary-state work and
    the legacy polarisation energy exchanges are separate, uncovered terms.
    """
    from .fluxes import fill_transport_guards
    from .inventory import energy_densities
    from .transport import fv_xz

    ds = run.data
    nu = run.scalar_option("viscosity", 0) / (run.length_scale**2 * run.frequency)
    vort = fill_transport_guards(ds.Vort, run, False, False)
    viscosity = xr.ones_like(ds.phi) * nu
    viscosity[{"x": run.mxg - 1}] = -nu
    viscosity[{"x": -run.mxg}] = -nu
    current = fv_xz(viscosity, vort, ds)
    Phi = ds.phi + energy_densities(run).pressure_potential

    def gradients(f):
        xp, xm = neighbor(f, "x"), neighbor(f, "x", -1)
        zp, zm = neighbor(f, "z"), neighbor(f, "z", -1)
        gx = (
            4 * (xp - xm)
            + neighbor(xp, "z")
            + neighbor(xp, "z", -1)
            - neighbor(xm, "z")
            - neighbor(xm, "z", -1)
        ) / (12 * ds.dx)
        gz = (
            4 * (zp - zm)
            + neighbor(zp, "x")
            + neighbor(zp, "x", -1)
            - neighbor(zm, "x")
            - neighbor(zm, "x", -1)
        ) / (12 * ds.dz)
        return gx, gz

    wx, wz = gradients(vort)
    px, pz = gradients(Phi)
    heat = -nu * (ds.g11 * wx * px + ds.g33 * wz * pz)
    kinetic = -Phi * current.divergence
    return xr.Dataset(
        {
            "DivJviscous": current.divergence,
            "kinematic_viscous_heating": heat,
            "native_viscous_kinetic_work": kinetic,
            "viscous_heating_defect": heat + kinetic,
        }
    )
