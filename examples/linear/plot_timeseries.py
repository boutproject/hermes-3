
from boutdata import collect
import numpy as np
import matplotlib.pyplot as plt

paths = ["annulus-isothermal-d"]

var = "Ne"
norm = collect("Nnorm", path=paths[0])
label = r"Electron density [m$^{-3}$]"
inds = (2,8,0)

datasets = [collect(var, path=path, xind=inds[0], yind=inds[1], zind=inds[2]).squeeze() for path in paths]
times = [collect("t_array", path=path) for path in paths]

data = np.concatenate(datasets) * norm
time = 1e3 * np.concatenate(times) / collect("Omega_ci", path=paths[0]) # ms


plt.plot(time, data)
plt.xlabel("Time [ms]")
plt.ylabel(label)
plt.savefig(var + "_{}_{}_{}_timeseries.png".format(*inds))
plt.savefig(var + "_{}_{}_{}_timeseries.pdf".format(*inds))
plt.show()

