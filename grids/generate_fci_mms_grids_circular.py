import os

import zoidberg as zb
import numpy as np


class Axial_circle(zb.field.Slab):
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

    inner = zb.rzline.shaped_line(R0=R0, a=rho_1, n=nz)
    outer = zb.rzline.shaped_line(R0=R0, a=rho_2, n=nz)

    pol_grid = zb.poloidal_grid.grid_annulus(inner, outer, nx + 4, nz, show=False)
    pol_grids = [pol_grid for _ in range(ny)]

    grid = zb.grid.Grid(
        pol_grids, np.linspace(0.0, 2.0 * np.pi, ny, endpoint=False), Ly, yperiodic=True
    )

    def calc_divertortheta(x, z, x0=0.0):
        theta = np.array(np.arctan2(z, x - x0))
        theta[theta < 0.0] += 2.0 * np.pi

        return theta

    with zb.zoidberg.MapWriter(fn) as mw:
        mw.add_grid_field(grid, field)

        maps = zb.zoidberg.make_maps(grid, field)
        mw.add_maps(maps)

        # Add some additional things
        rho = np.empty_like(maps["forward_xt_prime"])
        theta = np.empty_like(maps["forward_xt_prime"])
        phi = np.empty_like(maps["forward_xt_prime"])

        R = maps["R"][:, 0, :]
        Z = maps["Z"][:, 0, :]
        thistheta = calc_divertortheta(R, Z, x0=R0)
        R_mag = np.sqrt((R - R0) ** 2 + Z**2)
        rho[...] = ((R_mag - rho_1) / (rho_2 - rho_1))[:, None, :]
        theta[...] = thistheta[:, None, :]
        phi[...] = np.linspace(0, 2 * np.pi, ny, endpoint=False)[None, :, None]
        mw.write_dict(dict(rho=rho, theta=theta, phi=phi))

        mw.add_dagp()
