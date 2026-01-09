#!/usr/bin/env python3

import xhermes

# Following two lines needed so all variables are shown when printing the Dataset
import xarray as xr
xr.set_options(display_max_rows=1000)

# Set better figure size
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = (16,8)

ds = xhermes.open(".")
print(ds)

# Get rid of size-1 radial and toroidal directions
ds = ds.squeeze()

# Make an animation
# Note: saving a gif can be slow. Comment out `save_as` argument to disable.
ds.bout.animate_list(
    [["Ne", "Ni"], ["Ve", "Vi"], ["Pe", "Pi"]],
    show=True,
    save_as="hermes_animation",
)
