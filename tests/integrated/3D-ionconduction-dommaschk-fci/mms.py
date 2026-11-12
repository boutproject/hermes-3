#!/usr/bin/env python3

# Create the Dommaschk FCI grids for the ion heat conduction MMS test, and add
# the manufactured solution and the corresponding source to the grid files.
#
# The manufactured solution is defined in cylindrical coordinates (R, phi, Z)
# and the source is calculated exactly with sympy from the symbolic magnetic
# field of zoidberg's DommaschkPotentials. Both are evaluated on the grid points
# and read by Hermes-3 through [input:grid_variables].

import argparse
import os
import uuid
from pathlib import Path

import numpy as np
import zoidberg
from sympy import cos, diff, exp, lambdify, log, pi, sin, sqrt
from zoidberg.fieldtracer import FieldTracer, trace_poincare

parser = argparse.ArgumentParser()
parser.add_argument(
    "-v", "--verbose", help="increase output verbosity", action="store_true"
)
parser.add_argument(
    "-f", "--force", help="create grid files, even if they exist", action="store_true"
)
args = parser.parse_args()

orgprint = print


def print(*pargs, **kvargs):
    if args.verbose:
        orgprint(*pargs, **kvargs)


# Grid resolutions (nx, ny, nz), nx includes 2 x 2 guard cells
resolutions = [(20, 4, 32), (36, 8, 64), (68, 16, 128)]

# Dommaschk field of examples/stellarator-3D/dommaschk-3D/create_grid.py
C = np.zeros((6, 5, 4))
C[5, 2, 1] = -1.489
C[5, 2, 2] = -1.489
R0 = 2.0
Btor = 1.5
symmetry = 5
yperiod = 2.0 * np.pi / symmetry

# Inner and outer flux surfaces at the outboard midplane. These are nested
# surfaces inside the island chain, so that field lines only leave the domain
# through the first and last row, which are clamped below.
a1 = R0 + 0.03
a2 = R0 + 0.06

# Normalisations, as in BOUT.inp
Nnorm = 1e18
Tnorm = 5.0
Bnorm = 1.0
AA = 1.0

# Constants from bout/constants.hxx
qe = 1.602176634e-19
Mp = 1.672621898e-27
c = 299792458
mu0 = 4.0e-7 * np.pi
e0 = 1 / (c * c * mu0)

# Normalisations as in hermes-3.cxx
Cs0 = np.sqrt(qe * Tnorm / Mp)
Omega_ci = qe * Bnorm / Mp
rho_s0 = Cs0 / Omega_ci

field = zoidberg.field.DommaschkPotentials(C, R_0=R0, B_0=Btor)
R, phi, Z = field.R, field.phi, field.Z

##############################
# Manufactured solution (normalised), periodic with the field period

n = 1 + 0.2 * sin(5 * phi) * cos(10 * Z)
T = 1 + 0.4 * sin(5 * phi + 1) + 0.2 * cos(15 * (R - R0))
p = n * T

##############################
# Magnetic field. Same as in DommaschkPotentials, but symbolic. The constant
# normalisation to Btor does not change the unit vector or ratios of B.

BR = R0 * diff(field.P, R)
Bphi = R0 / R * diff(field.P, phi)
BZ = R0 * diff(field.P, Z)
B = sqrt(BR**2 + Bphi**2 + BZ**2)


def Grad_par(f):
    """b.Grad(f) in cylindrical coordinates"""
    return (BR * diff(f, R) + Bphi / R * diff(f, phi) + BZ * diff(f, Z)) / B


def Div_par(f):
    """Div(b f), using Div(B) = 0"""
    return B * Grad_par(f / B)


def softFloor(value, minimum):
    # As in hermes_utils.hxx, assuming value > 0
    return value + minimum * exp(-value / minimum)


##############################
# Ion heat conduction, as in braginskii_collisions.cxx and
# braginskii_conduction.cxx for a single ion species


def ion_ion_frequency(den, temp):
    """Ion-ion collision frequency [1/s], den in m^-3 and temp in eV"""
    Tlim = softFloor(temp, 0.1)
    # softFloor(den, 1e10) is den in double precision for den >> 1e10. Applying
    # it symbolically results in 0 * inf when evaluating the source
    Nlim = den
    coulomb_log = 29.91 - log((2 * AA) / (2 * AA * Tlim) * sqrt(2 * Nlim / Tlim))
    mass = AA * Mp
    vsq = 2 * Tlim * qe / mass
    return (
        qe**4
        * Nlim
        * softFloor(coulomb_log, 1.0)
        * 2
        / (3 * (pi * 2 * vsq) ** 1.5 * (e0 * mass) ** 2)
    )


tau = 1 / softFloor(ion_ion_frequency(n * Nnorm, T * Tnorm) / Omega_ci, 1e-10)
kappa = 3.9 * p * tau / AA

# Steady state: source balances conduction. Lengths are in m, so multiply by
# rho_s0^2 to get normalised units
Sp = -(2.0 / 3) * Div_par(kappa * Grad_par(T)) * rho_s0**2

lambdify_args = (R, phi, Z)
mms_functions = {
    "mms_Nh": lambdify(lambdify_args, n, "numpy"),
    "mms_Ph_solution": lambdify(lambdify_args, p, "numpy"),
    "mms_Ph_source": lambdify(lambdify_args, Sp, "numpy"),
}


def check_field():
    """Check the symbolic field and parallel gradient against zoidberg"""
    Rs = np.array([2.03, 2.05, 1.98])
    Zs = np.array([0.01, -0.02, 0.03])
    phis = np.array([0.1, 0.5, 1.0])

    # Same direction as zoidberg's field
    ratios = [
        lambdify(lambdify_args, Bsym, "numpy")(Rs, phis, Zs) / zfunc(Rs, Zs, phis)
        for Bsym, zfunc in (
            (BR, field.Bxfunc),
            (Bphi, field.Byfunc),
            (BZ, field.Bzfunc),
        )
    ]
    ratios = np.array(ratios)
    assert np.allclose(ratios, ratios[1, 0], rtol=1e-10), ratios

    # Grad_par(T) along a traced field line, compared to a central difference
    gradT = lambdify(lambdify_args, Grad_par(T) * B * R / Bphi, "numpy")
    Tf = lambdify(lambdify_args, T, "numpy")
    dphi = 1e-4
    tracer = FieldTracer(field)

    def follow(R1, Z1, phi1, phi2):
        """(R, Z) of the field line through (R1, phi1, Z1) at phi2"""
        return tracer.follow_field_lines(R1, Z1, [phi1, phi2], rtol=1e-12)[1, 0]

    for R1, Z1, phi1 in zip(Rs, Zs, phis):
        Rp, Zp = follow(R1, Z1, phi1, phi1 + dphi)
        Rm, Zm = follow(R1, Z1, phi1, phi1 - dphi)
        fd = (Tf(Rp, phi1 + dphi, Zp) - Tf(Rm, phi1 - dphi, Zm)) / (2 * dphi)
        assert np.isclose(fd, gradT(R1, phi1, Z1), rtol=1e-6), (fd, gradT(R1, phi1, Z1))
    print("Field checks passed")


def trace_surfaces(ny):
    """Trace the inner and outer flux surfaces on ny poloidal planes"""
    y_slices = np.linspace(0.0, yperiod, ny, endpoint=False)
    surfaces = []
    for a in (a1, a2):
        rzcoord, _ = trace_poincare(
            field, a, 0.0, yperiod, y_slices=y_slices, revs=500, nplot=1, nover=20
        )
        surfaces.append(rzcoord)
    return surfaces


def create_grid(nx, ny, nz):
    filename = Path(f"dommaschk_{nx}_{ny}_{nz}.fci.grid.nc")

    if filename.exists() and not args.force:
        print(filename, " exists")
        return

    # Note: The surfaces are traced for every grid, as subsampling the planes of
    # a finer grid results in a few field lines leaving the domain
    surfaces = trace_surfaces(ny)
    ycoords = np.linspace(0.0, yperiod, ny, endpoint=False)

    pol_grids = []
    for i in range(ny):
        inner, outer = [
            zoidberg.rzline.line_from_points(
                rz[:, i, 0, 0], rz[:, i, 0, 1], spline_order=3
            )
            for rz in surfaces
        ]
        pol_grids.append(
            zoidberg.poloidal_grid.grid_elliptic(
                inner,
                outer,
                nx,
                nz,
                restrict_size=2560,
                align=0,
                inner_ort=1,
                inner_maxmode=4,
                nx_inner=2,
                nx_outer=2,
            )
        )

    grid = zoidberg.grid.Grid(pol_grids, ycoords, yperiod, yperiodic=True)
    maps = zoidberg.make_maps(grid, field, nslice=1, quiet=not args.verbose)

    # The traced surfaces are slightly inaccurate, so field lines starting on the
    # first and last row can leave the domain. Bend them back into the domain,
    # as in examples/stellarator-3D/dommaschk-3D/create_grid.py
    for name in ("forward_xt_prime", "backward_xt_prime"):
        maps[name][2, :, :] = 2.0
        maps[name][-3, :, :] = nx - 3
        # No field line is allowed to leave the domain: the conduction operator
        # does not have boundary fluxes, which would be inconsistent with MMS
        interior = maps[name][2:-2]
        assert np.all((interior >= 2) & (interior <= nx - 3)), (
            f"{name}: field lines leave the domain in {filename}"
        )

    # Manufactured solution and source on the grid points
    R3d = maps["R"]
    Z3d = maps["Z"]
    phi3d = np.broadcast_to(ycoords[None, :, None], R3d.shape)
    mms_data = {name: func(R3d, phi3d, Z3d) for name, func in mms_functions.items()}
    assert np.min(mms_data["mms_Nh"]) * Nnorm > 1e15, "Density floor is active"

    tmpname = str(uuid.uuid4())
    with zoidberg.zoidberg.MapWriter(tmpname) as mw:
        mw.add_grid_field(grid, field)
        mw.add_maps(maps)
        mw.add_dagp()
        mw.write_dict(mms_data)
    os.replace(tmpname, filename)
    print(f"Created {filename}")


if __name__ == "__main__":
    check_field()

    for nx, ny, nz in resolutions:
        create_grid(nx, ny, nz)
