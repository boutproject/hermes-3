# Compare radial density profiles


import matplotlib.pyplot as plt
import numpy as np
from boutdata import collect

datasets = [
    ("D 64x64", "annulus-isothermal-d/", "b--"),
    ("D 128x128", "annulus-isothermal-d-2/", "b-"),
    ("He", "annulus-isothermal-he/", "r--"),
]

yind = 8

for label, path, style in datasets:
    n = collect("Ne", path=path, yind=yind).squeeze()
    nnorm = collect("Nnorm", path=path)
    nav = np.mean(n, axis=(0, 2)) * nnorm

    g_33 = collect("g_33", path=path, yind=yind).squeeze()
    radius = np.sqrt(g_33) * collect("rho_s", path=path)

    plt.plot(radius, nav, style, label=label)

plt.xlabel("Radius [m]")
plt.ylabel(r"Average plasma density [m$^{-3}$]")
plt.legend()

plt.savefig("profile-comparison.png")
plt.savefig("profile-comparison.pdf")

plt.show()
