import os

import zoidberg
import numpy as np
from zoidberg.field import Slab
import xarray as xr


class Axial_circle(Slab):
    def __init__(self, rho_1, rho_2, Lz, By):
        self.rho_1 = float(rho_1)
        self.rho_2 = float(rho_2)
        self.Lz = float(Lz)
        self.By = float(By)

    def Bxfunc(self, x, z, phi):
        return np.zeros(x.shape)

    def Bzfunc(self, x, z, phi):
        return np.zeros(x.shape)

    def Byfunc(self, x, z, phi):
        return np.full(x.shape, self.By)


todos = []

for scale in 1, 2, 4:
    nx = scale * 16
    ny = 4
    nz = scale * 64
    rho_1 = 0.2
    rho_2 = 0.4
    folder = "mms_circle_fci_xz"
    fn = f"{folder}/Axial_circular_{nx}_{nz}_{rho_1}_{rho_2}.fci.grid.nc"
    todos.append((nx, ny, nz, rho_1, rho_2, fn))

if not os.path.exists(folder):
    os.mkdir(folder)


for nx, ny, nz, rho_1, rho_2, fn in todos:
    if os.path.exists(fn):
        continue

    R0 = 1.5
    B0 = 2.5
    Ly = 2.0 * np.pi
    field = Axial_circle(rho_1, rho_2, 2.0 * np.pi / R0, B0)

    dx = (rho_2 - rho_1) / nx

    R1 = rho_1 - 1.5 * dx
    R2 = rho_2 + 1.5 * dx

    inner = zoidberg.rzline.shaped_line(R0=R0, a=rho_1, n=nz)
    outer = zoidberg.rzline.shaped_line(R0=R0, a=rho_2, n=nz)

    pol_grid = zoidberg.poloidal_grid.grid_annulus(inner, outer, nx + 4, nz)
    pol_grids = []
    for i in range(ny):
        pol_grids.append(pol_grid)

    grid = zoidberg.grid.Grid(
        pol_grids, np.linspace(0.0, 2.0 * np.pi, ny, endpoint=False), Ly, yperiodic=True
    )
    maps = zoidberg.zoidberg.make_maps(grid, field)

    zoidberg.zoidberg.write_maps(grid, field, maps, gridfile=fn, metric2d=False)

    gf = xr.open_dataset(fn)
    cp = gf.copy(deep=True)

    cp["rho"] = (cp["forward_xt_prime"].dims, np.zeros(cp["forward_xt_prime"].shape))
    cp["theta"] = (cp["forward_xt_prime"].dims, np.zeros(cp["forward_xt_prime"].shape))
    cp["phi"] = (cp["forward_xt_prime"].dims, np.zeros(cp["forward_xt_prime"].shape))

    def calc_divertortheta(x, z, x0=0.0):
        theta = np.array(np.arctan2(z, x - x0))
        theta[theta < 0.0] += 2.0 * np.pi

        return theta

    for i in range(cp["theta"].shape[0]):
        for k in range(cp["theta"].shape[2]):
            R = cp.R[i, 0, k].copy(deep=True).values
            Z = cp.Z[i, 0, k].copy(deep=True).values
            thistheta = calc_divertortheta(R, Z, x0=R0)
            R_mag = np.sqrt((R - R0) ** 2 + Z**2)
            rho = (R_mag - rho_1) / (rho_2 - rho_1)
            cp["theta"][i, :, k] = thistheta
            cp["rho"][i, :, k] = rho

    for i in range(cp["theta"].shape[1]):
        cp["phi"][:, i, :] = 2.0 * np.pi * i / ny

    cp.to_netcdf(fn)
