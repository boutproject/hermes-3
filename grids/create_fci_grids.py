import zoidberg
from zoidberg.field import Slab
import numpy as np

import os
import sys


force = "-f" in sys.argv or "--force" in sys.argv


class ThisField(Slab):
    def __init__(self, By, Byprime, xcentre, Bz=0.0):
        self.By = By
        self.Byprime = Byprime
        self.xcentre = xcentre
        self.Bz = Bz

    def Bxfunc(self, x, z, phi):
        return np.full(x.shape, 0.0)

    def Bzfunc(self, x, z, phi):
        return np.full(x.shape, self.Bz)

    def Byfunc(self, x, z, phi):
        return np.full(x.shape, self.By) + self.Byprime * np.sin(
            (x - self.xcentre) * 2.0 * np.pi
        )


class BlobField(Slab):
    def __init__(self, By, Byprime, xcentre, Bz=0.0):
        self.By = By
        self.Byprime = Byprime
        self.xcentre = xcentre
        self.Bz = Bz

    def Bxfunc(self, x, z, phi):
        return np.full(x.shape, 0.0)

    def Bzfunc(self, x, z, phi):
        return np.full(x.shape, self.Bz)

    def Byfunc(self, x, z, phi):
        return np.full(x.shape, self.By) + self.Byprime * (x - self.xcentre)


def create_directory(path):
    """Creates a directory if it does not already exist."""
    os.makedirs(path, exist_ok=True)
    print(f"Directory '{path}' is ready.")


def create_grid(folder, nx, ny, nz, BC=False, inp_Ly=None):

    filename = folder + f"MMS_straight_slab_{nx}_{ny}_{nz}_{BC}.fci.grid.nc"

    if os.path.exists(filename) and not force:
        print(filename, " exists")
        return 0

    rshift = 2.0
    magnetic_field = ThisField(By=1.0, xcentre=rshift, Byprime=-0.1)
    Lx = 1.0

    if BC:
        if inp_Ly is None:
            Ly = 100.0
        else:
            Ly = inp_Ly
    else:
        Ly = 2.0 * np.pi
    Lz = 1.0

    poloidal_grids = []
    for i in range(ny):
        tmp = zoidberg.poloidal_grid.RectangularPoloidalGrid(
            nx, nz, Lx, Lz, Rcentre=rshift
        )
        poloidal_grids.append(tmp)

    rectangle = zoidberg.grid.Grid(
        poloidal_grids, np.linspace(0, Ly, ny, endpoint=False), Ly, yperiodic=True
    )

    maps = zoidberg.make_maps(rectangle, magnetic_field, nslice=1)

    if BC:
        maps["forward_xt_prime"][:, -1, :] = nx - 1
        maps["backward_xt_prime"][:, 0, :] = nx - 1

    with zoidberg.zoidberg.MapWriter(filename) as mw:
        mw.add_grid_field(rectangle, magnetic_field)
        mw.add_maps(maps)
        mw.add_dagp()


script_dir = os.path.dirname(os.path.abspath(__file__))
print(script_dir)
folder_xz = script_dir + "/MMS_slab_xz/"
folder_y = script_dir + "/MMS_slab_y/"
folder_BC = script_dir + "/slab_with_BC/"


create_directory(folder_xz)
create_directory(folder_y)
create_directory(folder_BC)


create_grid(folder_y, 6, 16, 4)
create_grid(folder_y, 6, 32, 4)
create_grid(folder_y, 6, 64, 4)

create_grid(folder_xz, 36, 4, 32)
create_grid(folder_xz, 68, 4, 64)
create_grid(folder_xz, 132, 4, 128)

create_grid(folder_BC, 8, 128, 4, BC=True)
create_grid(folder_BC, 8, 512, 4, BC=True)
create_grid(folder_BC, 8, 256, 4, BC=True)


def create_blob_grid(nx, ny, nz, filename):

    if os.path.exists(filename) and not force:
        print(filename, " exists")
        return 0

    rshift = 1.5
    magnetic_field = BlobField(
        By=1.0 / rshift, xcentre=rshift, Byprime=-1.0 / (rshift**2)
    )
    Lx = 0.05
    Ly = 1.0
    Lz = 0.05

    poloidal_grids = []
    for i in range(ny):
        tmp = zoidberg.poloidal_grid.RectangularPoloidalGrid(
            nx, nz, Lx, Lz, Rcentre=rshift
        )
        poloidal_grids.append(tmp)

    rectangle = zoidberg.grid.Grid(
        poloidal_grids, np.linspace(0, Ly, ny, endpoint=False), Ly, yperiodic=True
    )

    maps = zoidberg.make_maps(rectangle, magnetic_field, nslice=1)

    with zoidberg.zoidberg.MapWriter(filename) as mw:
        mw.add_grid_field(rectangle, magnetic_field)
        mw.add_maps(maps)
        mw.add_dagp()


folder_blob = script_dir + "/slab_for_blob/"

create_directory(folder_blob)
create_blob_grid(132, 2, 128, folder_blob + "slab_blob_132_2_128.grid.fci.nc")
create_blob_grid(260, 2, 256, folder_blob + "slab_blob_260_2_256.grid.fci.nc")


print("Finished creating all the grids")
