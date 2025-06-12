#!/usr/bin/env python
#
# Makes plots of density and potential axisymmetric and non-axisymmetric components
# Use by passing the grid file and the directory with the BOUT.dmp.* files
#
# e.g.
#
#   python3 makeplots.py tokamak.nc isothermal/
#

import argparse

parser = argparse.ArgumentParser(description="Make plots of tokamak simulation output")
parser.add_argument("gridfile", type=str, help="The grid file used")
parser.add_argument("datapath", type=str, help="The path to the data (BOUT.dmp files)")

args = parser.parse_args()
gridfile = args.gridfile
path = args.datapath

from boutdata import collect
import matplotlib

#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter

cmap_diverging = plt.cm.get_cmap("bwr")  # RdBu, bwr, seismic
cmap_sequential = plt.cm.get_cmap("magma") # magma

from boututils.datafile import DataFile
from boutdata.griddata import gridcontourf
import numpy as np

#########################
# Read normalisations, time

Nnorm = collect("Nnorm", path=path)
Tnorm = collect("Tnorm", path=path)
wci = collect("Omega_ci", path=path)

t_arr = collect("t_array", path=path)
time = t_arr / wci # seconds

#######################
# Read the grid file

grid = DataFile(gridfile)

ixseps1 = grid["ixseps1"]
jyseps1_2 = grid["jyseps1_2"]
jyseps2_2 = grid["jyseps2_2"]
Rxy = grid["Rxy"]
# Poloidal index of midplane
j_midplane = jyseps1_2 + np.argmax(Rxy[ixseps1, jyseps1_2:(jyseps2_2+1)])

########################
# Read data to plot

ne = collect("Ne", path=path) * Nnorm  # m^-3
phi = collect("phi", path=path) * Tnorm  # V

# fix phi outer boundaries
phi[:,0,:,:] = 0.0
phi[:,-1,:,:] = 0.0

#######################
# Split into toroidal average and fluctuations

def acdc(f):
    """Calculate DC (Z average) and AC (Z varying) components
    """
    from numpy import mean, copy
    ftxy = mean(f, axis=-1)
    result = copy(f)
    for z in range(f.shape[-1]):
        result[:,:,:,z] -= ftxy
    return result, ftxy

ne_ac, ne_dc = acdc(ne)
phi_ac, phi_dc = acdc(phi)

ne_rms = np.sqrt(np.mean(ne_ac**2, axis=-1))
phi_rms = np.sqrt(np.mean(phi_ac**2, axis=-1))

######################

fig, ax = plt.subplots(1, 2, figsize=(10, 8))

c = gridcontourf(grid, ne_dc[-1,:,:], ax=ax[0], show=False, separatrix=True)
fig.colorbar(c, ax=ax[0])
ax[0].set_title(r"Electron density [m$^{-3}$]")
#ax[0].axes.get_xaxis().set_ticks([1, 2, 3, 4])

c = gridcontourf(grid, phi_dc[-1,:,:], ax=ax[1], show=False, ylabel=None, symmetric=True, cmap = cmap_diverging, separatrix=True)
fig.colorbar(c, ax=ax[1])
ax[1].set_title(r"Potential $\phi$ [V]")
ax[1].axes.get_yaxis().set_ticks([])
#ax[1].axes.get_xaxis().set_ticks([1, 2, 3, 4])

#fig.tight_layout(pad=3.0)

plt.savefig("ne_phi_axi.png")
plt.show()
plt.close("all")

#######################
# RMS fluctuation amplitude

fig, ax = plt.subplots(1, 2, figsize=(10, 8))

mind = 10.0**np.floor(np.log10(np.amin(ne_rms[-1,2:-2,:])))
maxd = 10.0**np.ceil(np.log10(np.amax(ne_rms[-1,2:-2,:])))

c = gridcontourf(grid, ne_rms[-1,:,:], ax=ax[0], show=False, log=True,
                 mind=mind, maxd=maxd, cmap=cmap_sequential)
cbar = fig.colorbar(c, ax=ax[0])
cbar.set_ticks(10.0**np.arange(np.log10(mind),np.log10(maxd)+1))
ax[0].set_title(r"RMS density $\sqrt{\left<\tilde{n_e}^2\right>}$ [m$^{-3}$]")
#ax[0].axes.get_xaxis().set_ticks([1, 2, 3, 4])


mind = 10.0**np.floor(np.log10(np.amin(phi_rms[-1,2:-2,:])))
maxd = 10.0**np.ceil(np.log10(np.amax(phi_rms[-1,2:-2,:])))

c = gridcontourf(grid, np.clip(phi_rms[-1,:,:], 1e-5, None), ax=ax[1], show=False,
                 ylabel=None, log=True, cmap=cmap_sequential,
                 mind=mind, maxd=maxd)
cbar = fig.colorbar(c, ax=ax[1])
cbar.set_ticks(10.0**np.arange(np.log10(mind),np.log10(maxd)+1))
ax[1].set_title(r"RMS potential $\sqrt{\left<\tilde{\phi}^2\right>}$ [V]")
ax[1].axes.get_yaxis().set_ticks([])
#ax[1].axes.get_xaxis().set_ticks([1, 2, 3, 4])

plt.savefig("ne_phi_rms.png")
plt.show()
plt.close("all")

###############

for j, label, short in [(0, "inner target", "inner"),
                        (j_midplane, "outboard midplane", "midplane"),
                        (-1, "outer target", "outer")]:
    data = ne_ac[-1,:,j,:].T

    mind = np.amin(data)
    maxd = np.amax(data)
    if maxd < -mind:
        maxd = -mind
    else:
        mind = -maxd
    levels = np.linspace(mind, maxd, 51, endpoint=True)
    plt.contourf(data, levels=levels, cmap=cmap_diverging)
    plt.axvline(ixseps1-0.5, linestyle='--', color='k')
    plt.colorbar()
    plt.title(r"Density $\tilde{n_e}$ at " + label)
    plt.xlabel("Radial index X")
    plt.ylabel("Toroidal index Z")
    plt.savefig("ne_ac_{}.png".format(short))
    plt.show()

#################

for j, label, short in [(0, "inner target", "inner"),
                        (j_midplane, "outboard midplane", "midplane"),
                        (-1, "outer target", "outer")]:
    data = np.amax(ne_rms[:,:,j], axis=-1)

    plt.plot(time*1e6, data, label=label)
plt.ylabel(r"Density $\tilde{n_e}$ at " + label)
plt.xlabel(r"Time [$\mu$s]")
plt.yscale('log')
plt.legend()
plt.savefig("ne_rms_time.png".format(short))
plt.show()
