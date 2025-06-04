# Calculate spatial correlation


from boutdata import collect
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import pickle

from scipy.signal import welch

paths = ["annulus-isothermal-d"]

var = "Ne"
label = r"Electron density [m$^{-3}$]"
norm = collect("Nnorm", path=paths[0])
two_sided = False

xind = 20
yind = 8

#var = "Vhe+"
#label = "Parallel flow speed [m/s]"
#norm = collect("Cs0", path=paths[0]) / 4  # divide by AA
#two_sided = True


rho_s = collect("rho_s", path=paths[0])
g_33 = collect("g_33", path=paths[0], yind=yind).squeeze()

if var[0] == "V":
    nvar = "N" + var[1:]
    
    datasets = [collect("N" + var, path=path, xind=xind, yind=yind).squeeze() / 
            np.clip(collect(nvar, path=path, xind=xind, yind=yind).squeeze(), 1e-5, None)
                for path in paths]
else:
    datasets = [collect(var, path=path, xind=xind, yind=yind).squeeze() for path in paths]

times = [collect("t_array", path=path) for path in paths]

data = np.concatenate(datasets) * norm
time = 1e3 * np.concatenate(times) / collect("Omega_ci", path=paths[0]) # ms

nt, nz = data.shape
print("Size: {}", data.shape)

print("Time range: {} - {}ms. Time samples: {}".format(time[0], time[-1], nt))

dt = time[1] - time[0]

f, psd = welch(data[:,0], fs=1./dt, nperseg=512) 
print(psd.shape)

for z in range(1, nz):
    _, p = welch(data[:,z], fs=1./dt, nperseg=512)
    psd += p
psd /= nz

import pickle
with open("psd.pkl", "wb") as pkl:
    pickle.dump(f, pkl)
    pickle.dump(psd, pkl)

plt.plot(f[1:], psd[1:])
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Frequency [kHz]')
plt.ylabel('Power Spectral Density')
plt.ylim([1e30, 1e38])

plt.savefig("psd-loglog.png")
plt.savefig("psd-loglog.pdf")

plt.show()

plt.plot(f[1:], psd[1:])
plt.yscale('log')
plt.xlim([0,20])
plt.xlabel('Frequency [kHz]')
plt.ylabel('Power Spectral Density')
plt.ylim([1e30, 1e38])

plt.savefig("psd-loglin.png")
plt.savefig("psd-loglin.pdf")

plt.show()

