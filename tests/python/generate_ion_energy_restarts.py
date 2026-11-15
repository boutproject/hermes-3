"""Generate smooth single-rank X-Z restart files with constant potential shifts."""

from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def generate(folder, shift, n=8):
    """Write identical physical evolved cells; encode the gauge in wall storage.

    Vorticity guards store the fixed wall potential. Potential guards reflect
    about that trace, while physical vorticity, density, pressure and energy
    remain independent of the shift. No preparation simulation is needed.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    x = (np.arange(n + 4) - 1.5) / n
    z = 2 * np.pi * np.arange(n) / n
    X, Z = np.meshgrid(x, z, indexing="ij")
    fields = {
        "Vort": 0.001 * (1 + 0.2 * np.cos(2 * np.pi * X)) * np.sin(Z),
        "Nh+": 1 + 0.05 * np.sin(2 * np.pi * X) * np.cos(Z),
        "Pe": 1 + 0.08 * np.sin(np.pi * X) * np.sin(Z),
        "Eion": 10 + 0.05 * np.cos(2 * np.pi * X) * np.cos(Z),
        "phi": 0.04 * np.sin(np.pi * X) * np.sin(Z) + shift,
    }
    for name, field in fields.items():
        if name == "Vort":
            field[:2] = field[-2:] = shift
        elif name == "phi":
            field[:2] = 2 * shift - field[2]
            field[-2:] = 2 * shift - field[-3]
        else:
            field[:2] = field[2]
            field[-2:] = field[-3]
    path = folder / "BOUT.restart.0.nc"
    with Dataset(path, "w", format="NETCDF4_CLASSIC") as dump:
        for name, size in (("x", n + 4), ("y", 1), ("z", n)):
            dump.createDimension(name, size)
        for name, field in fields.items():
            variable = dump.createVariable(name, "f8", ("x", "y", "z"))
            variable[:] = field[:, None, :]
            variable.setncatts(
                {
                    "cell_location": "CELL_CENTRE",
                    "direction_y": "Standard",
                    "direction_z": "Standard",
                }
            )
        scalars = {
            "tt": 0.0,
            "hist_hi": 0,
            "MXG": 2,
            "MYG": 0,
            "MZG": 0,
            "NXPE": 1,
            "NYPE": 1,
            "NZPE": 1,
            "MYPE": 0,
            "MXSUB": n,
            "MYSUB": 1,
            "MZSUB": n,
            "MZ": n,
            "PE_XIND": 0,
            "PE_YIND": 0,
            "nx": n + 4,
            "ny": 1,
            "nz": n,
            "ion_energy_boundary_state_version": 1,
        }
        for name, value in scalars.items():
            variable = dump.createVariable(
                name, "i4" if isinstance(value, int) else "f8"
            )
            variable.assignValue(value)
    return path


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    for shift in (-1.0, 0.0, 1.0):
        generate(args.folder / f"shift{shift:g}", shift)
