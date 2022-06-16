from boutdata import collect
import numpy as np
import matplotlib.pyplot as plt

path = "annulus-isothermal-d"
var = "Ne"
norm = collect("Nnorm", path=path)
label= r"Electron density [m$^{-3}$]"
yind = 8
tind = -1

g_33 = collect("g_33", path=path, yind=yind).squeeze()
radius = np.sqrt(g_33) * collect("rho_s", path=path)
dz = collect("dz", path=path)[0,0]

data = collect(var, path=path, tind=tind, yind=yind).squeeze()

nz = data.shape[-1]
angle = np.arange(nz) * dz

r, theta = np.meshgrid(radius, angle, indexing='ij')

x = r * np.cos(theta)
y = r * np.sin(theta)

def close_period(arr):
    return np.concatenate([arr, arr[:,0].reshape((arr.shape[0], 1))], axis=1)

plt.contourf(close_period(x), close_period(y), close_period(data), 50)
plt.colorbar()
plt.xlabel("[m]")
plt.ylabel("[m]")
#plt.axes().set_aspect('equal')
plt.savefig(var + "_crossection.png")
plt.savefig(var + "_crossection.pdf")
plt.show()

