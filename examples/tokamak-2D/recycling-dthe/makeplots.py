#!/usr/bin/env python
#
# Makes plots of density and temperatures axisymmetric components
# Use by passing the grid file and the directory with the BOUT.dmp.* files
#
# e.g.
#
#   python3 makeplots.py tokamak.nc recycling-dthe
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
cs = collect("Cs0", path=path)

t_arr = collect("t_array", path=path)
time = t_arr / wci # seconds

#######################
# Read the grid file

grid = DataFile(gridfile)

ixseps1 = grid["ixseps1"]
jyseps1_2 = grid["jyseps1_2"]
jyseps2_2 = grid["jyseps2_2"]
Rxy = grid["Rxy"]
Zxy = grid["Zxy"]
# Poloidal index of midplane
j_midplane = jyseps1_2 + np.argmax(Rxy[ixseps1, jyseps1_2:(jyseps2_2+1)])

########################
# Read data to plot

datasets = [("Ne", Nnorm, r"Electron density [m$^{-3}$]", '-')
           ,("Te", Tnorm, r"Electron temperature [eV]", '-')
           ,("Td+", Tnorm, r"D+ ion temperature [eV]", '--')
           ,("Tt+", Tnorm, r"T+ ion temperature [eV]", '-.')
           ,("The+", Tnorm, r"He+ ion temperature [eV]", ':')
           ,("Nd", Nnorm, r"D atom density [m$^{-3}$]", '--')
           ,("Nt", Nnorm, r"T atom density [m$^{-3}$]", '-.')
           ,("Nhe", Nnorm, r"He atom density [m$^{-3}$]", ':')
           #,("Vd+", cs, r"Deuterium ion parallel flow [m/s]")
           ]

alldata = {}

for name, norm, _, _ in datasets:
  if name[0] == "T":
    p = collect("P" + name[1:], path=path, tind=-1).squeeze()
    n = collect("N" + name[1:], path=path, tind=-1).squeeze().clip(min=1e-5)
    alldata[name] = (p / n) * Tnorm
  elif name[0] == "V":
    nv = collect("N" + name, path=path, tind=-1).squeeze()
    n = collect("N" + name[1:], path=path, tind=-1).squeeze().clip(min=1e-5)
    alldata[name] = (nv / n) * Tnorm
  else:
    alldata[name] = collect(name, path=path, tind=-1).squeeze() * norm

nhe1 = collect("Nhe+", path=path,tind=-1).squeeze() * Nnorm # Helium ion density

datasets.append(("Che", 1.0, "Helium concentration (%)", ':'))
alldata["Che"] = 100. * nhe1 / alldata["Ne"]

######################

for name, _, label, _ in datasets:
  data = alldata[name]

  fig, ax = plt.subplots(1, 1, figsize=(10, 8))

  if name[0] == "V":
    c = gridcontourf(grid, data[:,:], ax=ax, show=False, symmetric=True, cmap = cmap_diverging, separatrix=True, remove_xguards=True)
  else:
    c = gridcontourf(grid, data[:,:], ax=ax, show=False, separatrix=True, remove_xguards=True)
  fig.colorbar(c, ax=ax)
  ax.set_title(label)
  #ax.axes.get_xaxis().set_ticks([1, 2, 3, 4])

  #fig.tight_layout(pad=3.0)

  plt.savefig(name + ".png")
  plt.savefig(name + ".pdf")
  plt.show()
  plt.close("all")

#####################
# Radial profiles

def replace_guards(arr1d):
    result = arr1d[1:-1].copy()
    result[0] = 0.5*(result[0] + result[1])
    result[-1] = 0.5*(result[-1] + result[-2])
    return result

for j, location in [(j_midplane, "midplane"), (0, "inner target"), (-1, "outer target")]:

  fig, ax = plt.subplots(1, 1, figsize=(10, 8))

  dR = Rxy[1:, j] - Rxy[:-1, j]
  dZ = Zxy[1:, j] - Zxy[:-1, j]
  dl = np.sqrt(dR**2 + dZ**2)
  distance = np.concatenate(([0], np.cumsum(dl))) * 1e2 # cm
  distance -= 0.5*(distance[ixseps1] + distance[ixseps1-1]) # From separatrix
  distance = replace_guards(distance)

  for name, _, label, linestyle in datasets:
    if name[0] != "N":
      continue
    data = replace_guards(alldata[name][:,j])
    ax.plot(distance, data, label=label, color='b', linestyle=linestyle)
  ax.legend(loc='lower left' if location == "midplane" else 'upper left')
  ax.set_ylim(bottom=0.0)

  # Plot temperatures
  ax2 = ax.twinx()
  for name, _, label, linestyle in datasets:
    if name[0] != "T":
      continue
    data = replace_guards(alldata[name][:,j])
    ax2.plot(distance, data, label=label, color='r', linestyle=linestyle)
  ax2.legend(loc='upper right')
  ax2.set_ylim(bottom = 0.0)

  ax.set_xlabel("Distance from separatrix [cm]")
  ax.axvline(0.0, linestyle='--', color='k')

  plt.savefig(location + ".png")
  plt.savefig(location + ".pdf")
  plt.show()
  plt.close("all")


