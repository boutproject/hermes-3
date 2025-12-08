#!/usr/bin/env python3

import xhermes

# Following two lines needed so all variables are shown when printing the Dataset
import xarray as xr
xr.set_options(display_max_rows=1000)

# Set better figure size
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = (16,8)

ds = xhermes.open(".", geometry="toroidal", gridfilepath="tokamak.nc")
print(ds)

# Get rid of size-1 toroidal direction
ds = ds.squeeze()

# Make some animations
# Note: saving a gif can be slow. Comment out `save_as` argument to disable.
ds.bout.animate_list(
    ["Ne", "Ve", "Pe"],
    poloidal_plot=True,
    ncols=3,
    show=True,
    save_as="electrons",
)
ds.bout.animate_list(
    ["Nd+", "NVd+", "Pd+", "Nd", "NVd", "Pd"],
    poloidal_plot=True,
    ncols=3,
    show=True,
    save_as="deuterium",
)
ds.bout.animate_list(
    ["Nt+", "NVt+", "Pt+", "Nt", "NVt", "Pt"],
    poloidal_plot=True,
    ncols=3,
    show=True,
    save_as="tritium",
)
ds.bout.animate_list(
    ["Nhe+", "NVhe+", "Phe+", "Nhe", "NVhe", "Phe"],
    poloidal_plot=True,
    ncols=3,
    show=True,
    save_as="helium",
)
