#!/usr/bin/env python

# Makes plots of the results

import matplotlib


import matplotlib.pyplot as plt
from boutdata import collect
from boututils.datafile import DataFile

import numpy as np
import csv

gridfilename = "../uedge.grd_Up_Ni_Tei_2d.nc"

def readcsv(csvfile):
    """
    Assumes a single header line, followed by columns of real numbers
    """
    reader = csv.reader(csvfile, delimiter=',')
    
    header = next(reader)
    data = []
    for row in reader:
        rowdata = [float(val) for val in row]
        data.append(rowdata)
    return header, np.array(data)

# Read UEDGE benchmark data
with open("ue_bmk.csv") as f:
    _, ue_data = readcsv(f)

# Change time base to ms
ue_data[:,0] *= 1e3

# Read grid
with DataFile(gridfilename) as d:
    Rxy = d["Rxy"]
    Zxy = d["Zxy"]

# Collect data

time = collect("t_array")
nt = len(time)  # Make sure that all fields have the same size

n = collect("Ne", tind=[0,nt-1])
nvi = collect("NVi", tind=[0,nt-1])
ti = collect("Ti", tind=[0,nt-1])
te = collect("Te", tind=[0,nt-1])


vi = nvi / n

# Read normalisation
Tnorm = collect("Tnorm")
Nnorm = collect("Nnorm")
Cs0 = collect("Cs0")
wci = collect("Omega_ci")

# Time into ms
time *= 1e3 / wci

# (x,y) indices where time traces should be plotted
#locations = [ (9,39), (19,39), (39,39) ]
#locations = [ (9,40), (19,40), (39,40) ]
locations = [ (10,40), (20,40), (40,40) ]
#locations = [(11,40), (21,40), (41,40) ]

ue_color = '#E0E0E0'

##########################################################

plt.figure()

# Density

plt.subplot(221)

plt.plot(ue_data[:,0], ue_data[:,1], color=ue_color)
plt.plot(ue_data[:,0], ue_data[:,2], color=ue_color)
plt.plot(ue_data[:,0], ue_data[:,3], color=ue_color)

for xi,yi in locations:
    plt.plot(time, n[:,xi,yi,0]*Nnorm)
#plt.xlabel("Time [ms]")
plt.title(r"Density [$m^{-3}$]")


# Vi

plt.subplot(222)

plt.plot(ue_data[:,0], ue_data[:,10]*1e-3, color=ue_color)
plt.plot(ue_data[:,0], ue_data[:,11]*1e-3, color=ue_color)
plt.plot(ue_data[:,0], ue_data[:,12]*1e-3, color=ue_color)

for xi,yi in locations:
    plt.plot(time, vi[:,xi,yi,0]*Cs0*1e-3)
#plt.xlabel("Time [ms]")
plt.title("Ion velocity [km/s]")

# Ti

plt.subplot(223)

plt.plot(ue_data[:,0], ue_data[:,7], color=ue_color)
plt.plot(ue_data[:,0], ue_data[:,8], color=ue_color)
plt.plot(ue_data[:,0], ue_data[:,9], color=ue_color)

for xi,yi in locations:
    plt.plot(time, ti[:,xi,yi,0]*Tnorm)
plt.xlabel("Time [ms]")
plt.title("Ion temperature [eV]")

# Te

plt.subplot(224)

plt.plot(ue_data[:,0], ue_data[:,4], color=ue_color)
plt.plot(ue_data[:,0], ue_data[:,5], color=ue_color)
plt.plot(ue_data[:,0], ue_data[:,6], color=ue_color)

for xi,yi in locations:
    plt.plot(time, te[:,xi,yi,0]*Tnorm)
plt.xlabel("Time [ms]")
plt.title("Electron temperature [eV]")


plt.savefig("hermes-history.pdf")

##########################################################

from boutdata.griddata import gridcontourf
from boututils.datafile import DataFile

matplotlib.rcParams.update({'font.size': 8})

grid = DataFile(gridfilename)

fig, axarr = plt.subplots(2,2)

# Density

axis = axarr[0,0]
c = gridcontourf(grid, n[-1,:,:,0]*Nnorm, show=False, ax=axis,cmap=plt.cm.get_cmap("Reds"), xlabel=None)
fig.colorbar(c, ax=axis)
axis.set_title(r"Density [$m^{-3}$]")

# Vi
axis = axarr[0,1]
c = gridcontourf(grid, vi[-1,:,:,0]*Cs0*1e-3, show=False, ax=axis, cmap=plt.cm.get_cmap("seismic"), symmetric=True, xlabel=None, ylabel=None)
fig.colorbar(c, ax=axis)
axis.set_title("Ion velocity [km/s]")

# Ti
axis = axarr[1,0]
c = gridcontourf(grid, ti[-1,:,:,0]*Tnorm, show=False, ax=axis)
fig.colorbar(c, ax=axis)
axis.set_title("Ion temperature [eV]")

# Te
axis = axarr[1,1]
c = gridcontourf(grid, te[-1,:,:,0]*Tnorm, show=False, ax=axis, ylabel=None)
fig.colorbar(c, ax=axis)
axis.set_title("Electron temperature [eV]")

plt.savefig("hermes-results.pdf")

plt.show()

matplotlib.rcParams.update({'font.size': 12})

vi = collect("vi", yguards=True, tind=-1, xind=25)
vi = vi[-1,0,1:-1,0] * Cs0

n = collect("Ne", yguards=True, tind=-1, xind=25)
n = n[-1,0,1:-1,0] * Nnorm

ny = len(vi)
yind = np.arange(ny)-1

fig, ax1 = plt.subplots()

ax1.plot(yind, vi*1e-3, 'k')
ax1.set_xlabel("Y index")
ax1.set_ylabel(r"Parallel velocity $v_{||i}$ [km/s]")
plt.axvline(-0.5, ls='--')
plt.axvline(ny-2.5, ls='--')

plt.title(r"Time = $%.1f$ ms, $x=25$" % (time[-1]))

ax1.set_xlim([-1,64])

ax2 = ax1.twinx()

ax2.plot(yind, n, 'r')
ax2.set_ylabel(r"Density [m$^{-3}$]")
ax2.tick_params('y', colors='r')
fig.tight_layout()

plt.savefig("hermes-parallel.pdf")

plt.show()
