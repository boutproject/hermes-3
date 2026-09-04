#!/usr/bin/env python3

# Visualise results from the 1D-low-n-srcs test:
#  1. An animation of the density evolution
#  2. A plot of the initial densities.

from __future__ import division
from __future__ import print_function

import os.path
import pathlib
import xhermes

import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter

this_dir = pathlib.Path(__file__).parent
data_dir = this_dir / "data"
hermes_exec = this_dir.parent.parent.parent / "hermes-3"
data_path = data_dir / "BOUT.dmp.0.nc"
scale_anim_size = 1.2

if not data_path.is_file():
    print(f"No data to plot at [{data_path}], bailing out")
    exit(1)

# Read data
ds = xhermes.open(data_dir, unnormalise=False)

# Species densities
species_list = ["d", "d+", "t", "t+"]
dens_group_indices = [[0, 1], [2, 3]]
dens_var_groups = [
    [f"N{species_list[i]}" for i in group] for group in dens_group_indices
]
dens_subplot_titles = ["Deuterium atoms,ions", "Tritium atoms,ions"]

# Source diagnostics
src_list = [
    "Sd_low-n_recombination",
    "Sd+_low-n_ionization",
    "Sd+_low-n_external_src",
    "St_low-n_recombination",
    "St+_low-n_ionization",
    "St+_low-n_external_src",
]
src_group_indices = [[0, 1, 2], [3, 4, 5]]
src_var_groups = [[f"{src_list[i]}" for i in group] for group in src_group_indices]
anim_vars = [[f"N{species_list[i]}" for i in group] for group in dens_group_indices]
anim_vars.extend(src_var_groups)
anim_subplot_titles = list(dens_subplot_titles)
anim_subplot_titles.extend(["Deuterium sources", "Tritium sources"])
anim_logscales = [True, True, False, False]

# Extract density thresholds
options = ds.options
thresh_fac = options["low_n_sources"]["low_n_floor_fac"]
floors = {}
thresholds = {}
for species in species_list:
    if "density_floor" not in options[species]:
        print(f"WARNING: Didn't find density floor for {species}; guessing 1e-7")
        floor = 1e-7
    else:
        floor = options[species].as_dict().get("density_floor", 1e-7)
    floors[species] = floor
    thresholds[species] = floor * thresh_fac


print("Generating ICs plot")
ICs_output_path = os.path.realpath("low-n_sources_initial")
ds_first = ds.isel(t=0)

fig, axes = plt.subplots(len(dens_var_groups), 1, sharex=True)
for ax, title, group, species_group in zip(
    axes, dens_subplot_titles, dens_var_groups, dens_group_indices
):
    for var in group:
        # animate_list() mutates density_vars in place, replacing names with DataArrays
        var_name = var if isinstance(var, str) else var.name
        ds_first[var_name].plot(ax=ax, label=var_name)
    ax.set_yscale("log")
    ax.set_ylim(bottom=min(floors[species_list[i]] / 2 for i in species_group))
    for i in species_group:
        species = species_list[i]
        color = ax.get_lines()[species_group.index(i)].get_color()
        linestyle = ":" if species.endswith("+") else "--"
        ax.axhline(
            thresholds[species],
            color=color,
            linestyle=linestyle,
            label=f"N{species}_th",
        )
    ax.set_ylabel("number density [normalised]")
    ax.set_title(title)
    ax.legend()
fig.tight_layout()
fig.savefig(ICs_output_path + ".png")
print(f"Wrote output to {ICs_output_path}.png")

print("Generating animation")
anim_output_path = os.path.realpath("low-n_sources_animation")
dens_vmin_values = [floors[species] for species in species_list]
anim_vmins = [min(dens_vmin_values[i] for i in group) for group in dens_group_indices]
anim_vmins.extend([0] * len(src_var_groups))  # Sources can go to zero
animation = ds.bout.animate_list(
    anim_vars,
    vmin=anim_vmins,
    titles=anim_subplot_titles,
    logscale=anim_logscales,
    show=False,
    save_as=None,
    controls=None,
)

# Increase the size of the exported animation figure
anim_fig = animation.blocks[0].ax.figure
fig_width, fig_height = anim_fig.get_size_inches()
anim_fig.set_size_inches(
    fig_width * scale_anim_size, fig_height * scale_anim_size, forward=True
)

# Remove superfluous color bar axes
plot_axes = {block.ax for block in animation.blocks}
for axis in list(animation.blocks[0].ax.figure.axes):
    if axis in plot_axes:
        axis.set_xlabel("Parallel connection length [normalised]")
    else:
        axis.remove()

# animation.blocks has one entry per line, in the order the variables were
# passed to animate_list, i.e. all density panel lines followed by all source panel lines
n_dens_blocks = sum(len(group) for group in dens_var_groups)

# Add horizontal lines to indicate density thresholds, colored consistently with the corresponding density variable
for block, species in zip(animation.blocks[:n_dens_blocks], species_list):
    color = block.line.get_color()
    linestyle = ":" if species.endswith("+") else "--"
    block.ax.axhline(
        thresholds[species],
        color=color,
        linestyle=linestyle,
        label=f"N{species}_th",
    )
    block.ax.legend(loc="lower left", fontsize=8)
    current_ylim = block.ax.get_ylim()
    block.ax.set_ylabel("[normalised] density")
    block.ax.set_ylim(
        bottom=min(floors[species] / 2, current_ylim[0]),
        top=max(thresholds[species] * 2, current_ylim[1]),
    )

for block in animation.blocks[n_dens_blocks:]:
    block.ax.set_ylabel("[normalised] source")
    block.ax.legend(fontsize=8)

animation.save(anim_output_path + ".gif", writer=PillowWriter(fps=10))
print(f"Wrote output to {anim_output_path}.gif")
