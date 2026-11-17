#!/usr/bin/env python3
"""Create the unstructured FCI grids used by the Hermes-3 tests.

Each grid is the unstructured counterpart of a structured grid made by
``grids/create_fci_grids.py``: same magnetic field, same domain, and
poloidal planes that are triangle meshes stored "nodes along x" (nz = 1).
See ``zoidberg.unstructured`` for the grid-file contract.

    python create_unstructured_grids.py            # Hermes-3 test grids
    python create_unstructured_grids.py --bout-mms # also BOUT-dev MMS grids

Existing files are kept unless ``-f`` is given. Every grid uses a fixed
seed, so regenerating gives identical data (only the random grid_id
attribute differs).
"""

import argparse
import os
import uuid
from pathlib import Path

import numpy as np
from zoidberg.field import Slab
from zoidberg.unstructured import (
    UnstructuredPoloidalGrid,
    make_test_grid,
    make_unstructured_grid,
)


class ThisField(Slab):
    """Same field as ``ThisField`` in ``grids/create_fci_grids.py``.

    Copied here because that script parses its arguments at import time.
    """

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


# Structured grids: RectangularPoloidalGrid(nx, nz, Lx=1, Lz=1, Rcentre=2)
RSHIFT = 2.0
LX = 1.0
LZ = 1.0


def field():
    return ThisField(By=1.0, xcentre=RSHIFT, Byprime=-0.1)


def z0(nz):
    """Z of the first point of the structured grid with nz points, so that
    BOUT++'s z = 2 pi (Z - z0) / Lz matches the structured grid."""
    return -0.5 * LZ + 0.5 * LZ / nz


def slab_y_poloidal():
    """Counterpart of MMS_slab_y (nx = 6 with MXG = 2, nz = 4): 6 x 8 nodes
    on R in [1.5, 2.5], periodic in Z. The field is straight, so every
    parallel map is an injection."""
    return UnstructuredPoloidalGrid.slab(
        6,
        8,
        LX,
        LZ,
        R0=RSHIFT - 0.5 * LX,
        Z0=z0(8),
        jitter=0.2,
        jitter_dirs="z",
        periodic_R=False,
        periodic_Z=True,
        seed=100,
    )


def slab_xz_poloidal(n):
    """Counterpart of MMS_slab_xz (n x n interior points): n x n nodes,
    doubly periodic. The MMS solutions are periodic in x and z, and
    ThisField is periodic in R with period Lx = 1."""
    return UnstructuredPoloidalGrid.slab(
        n,
        n,
        LX,
        LZ,
        R0=RSHIFT - 0.5 * LX,
        Z0=z0(n),
        jitter=0.2,
        jitter_dirs="rz",
        periodic_R=True,
        periodic_Z=True,
        seed=200 + n,
    )


def slab_xz_walls_poloidal(n):
    """Counterpart of MMS_slab_xz with its Dirichlet x walls: n + 1 node
    layers from R = 1.5 to 2.5 (nodes on both walls), periodic in Z. The wall
    nodes become Dirichlet nodes (dirichlet_walls=True)."""
    return UnstructuredPoloidalGrid.slab(
        n + 1,
        n,
        LX,
        LZ,
        R0=RSHIFT - 0.5 * LX,
        Z0=z0(n),
        jitter=0.2,
        jitter_dirs="rz",
        periodic_R=False,
        periodic_Z=True,
        seed=300 + n,
        r_on_walls=True,
    )


def write(filename, force, build):
    if filename.exists() and not force:
        print(filename, "exists")
        return
    tmpname = filename.parent / f".{uuid.uuid4()}.nc"
    build(str(tmpname))
    os.replace(tmpname, filename)
    print("Created", filename)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-f", "--force", action="store_true", help="overwrite grids")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument(
        "--bout-mms",
        action="store_true",
        help="also create the BOUT-dev MMS-spatial-fci-unstructured grids",
    )
    parser.add_argument(
        "--only-bout-mms", action="store_true", help="only create the BOUT-dev grids"
    )
    args = parser.parse_args()

    here = Path(__file__).resolve().parent

    if not args.only_bout_mms:
        folder = here / "US_slab_y"
        folder.mkdir(exist_ok=True)
        for ny in (16, 32, 64):
            write(
                folder / f"US_slab_y_{ny}.fci.grid.nc",
                args.force,
                lambda fn, ny=ny: make_unstructured_grid(
                    fn, slab_y_poloidal(), field(), ny=ny, Ly=2.0 * np.pi
                ),
            )

        folder = here / "US_slab_xz"
        folder.mkdir(exist_ok=True)
        for n in (16, 32):
            write(
                folder / f"US_slab_xz_{n}.fci.grid.nc",
                args.force,
                lambda fn, n=n: make_unstructured_grid(
                    fn, slab_xz_poloidal(n), field(), ny=1, Ly=2.0 * np.pi
                ),
            )

        folder = here / "US_slab_xz_walls"
        folder.mkdir(exist_ok=True)
        for n in (16, 32):
            write(
                folder / f"US_slab_xz_walls_{n}.fci.grid.nc",
                args.force,
                lambda fn, n=n: make_unstructured_grid(
                    fn,
                    slab_xz_walls_poloidal(n),
                    field(),
                    ny=1,
                    Ly=2.0 * np.pi,
                    dirichlet_walls=True,
                ),
            )

        folder = here / "US_slab_with_BC"
        folder.mkdir(exist_ok=True)
        for ny in (256,):
            write(
                folder / f"US_slab_with_BC_{ny}.fci.grid.nc",
                args.force,
                lambda fn, ny=ny: make_unstructured_grid(
                    fn, slab_y_poloidal(), field(), ny=ny, Ly=100.0, walls_y=True
                ),
            )

    if args.bout_mms or args.only_bout_mms:
        folder = here / "bout_mms"
        folder.mkdir(exist_ok=True)
        for n in (8, 16, 32, 64):
            write(
                folder / f"US_fci_mms_{n}.fci.nc",
                args.force,
                lambda fn, n=n: make_test_grid("bout_mms", n, filename=fn),
            )

    print("Finished creating unstructured grids")


if __name__ == "__main__":
    main()
