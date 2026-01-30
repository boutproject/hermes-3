import os
import sys

import numpy as np

import zoidberg


folder = "mms_slab_fci_xz/"
todo = []
for scale in 2, 4:
    nx = 16 * scale + 4  # 4 Guard cells
    nz = 16 * scale
    ny = 4
    todo.append((folder, nx, ny, nz))
folder = "mms_slab_fci_y/"
for scale in 1, 2, 4, 8:
    nx = 2 + 4  # 4 Guard cells
    ny = 16 * scale
    nz = 4
    todo.append((folder, nx, ny, nz))

force = "-f" in sys.argv or "--force" in sys.argv

for folder, nx, ny, nz in todo:
    filename = f"MMS_straight_slab_{nx}_{ny}_{nz}.fci.grid.nc"
    fn = folder + filename

    if os.path.exists(fn) and not force:
        continue

    rshift = 1.5
    magnetic_field = zoidberg.field.Slab(By=1.0, Bz=0.0, Bzprime=0.0, xcentre=rshift)
    Lx = 1.0
    Ly = 2.0 * np.pi
    Lz = 1.0

    poloidal_grids = [
        zoidberg.poloidal_grid.RectangularPoloidalGrid(nx, nz, Lx, Lz, Rcentre=rshift)
        for _ in range(ny)
    ]

    rectangle = zoidberg.grid.Grid(
        poloidal_grids, np.linspace(0, Ly, ny, endpoint=False), Ly, yperiodic=True
    )

    qq = rectangle.metric()
    maps = zoidberg.make_maps(rectangle, magnetic_field, nslice=2)

    with zoidberg.zoidberg.MapWriter(fn) as mw:
        mw.add_grid_field(rectangle, magnetic_field)
        mw.add_maps(maps)
        mw.add_dagp()
    print(fn, "done")
