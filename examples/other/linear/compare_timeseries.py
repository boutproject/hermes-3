
from boutdata import collect
import numpy as np
import matplotlib.pyplot as plt
import os

datasets = [
    ("D", "annulus-isothermal-d", ["."])
    ,("He", "annulus-isothermal-he", ["."])
           ]

var = "Ne"
ylabel = r"Electron density [m$^{-3}$]"
inds = (2,8,0)

for label, base, subdirs in datasets:
    paths = [os.path.join(base, subdir) for subdir in subdirs]
    
    norm = collect("Nnorm", path=paths[0])
    datas = [collect(var, path=path, xind=inds[0], yind=inds[1], zind=inds[2]).squeeze() for path in paths]
    times = [collect("t_array", path=path) for path in paths]

    data = np.concatenate(datas) * norm
    time = 1e3 * np.concatenate(times) / collect("Omega_ci", path=paths[0]) # ms


    plt.plot(time, data, label=label)

plt.xlabel("Time [ms]")
plt.ylabel(ylabel)
plt.legend()

plt.savefig(var + "_{}_{}_{}_comparison.png".format(*inds))
plt.savefig(var + "_{}_{}_{}_comparison.pdf".format(*inds))
plt.show()

