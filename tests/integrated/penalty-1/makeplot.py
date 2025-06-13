#!/usr/bin/env python3

from boutdata import collect
import numpy as np

import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 2)
axs = axs.flatten()

varlist = ["Nd+", "NVd+", "Pd+", "Pe"]

for path in ["sheath", "penalty", "penalty-e", "sheath-e"]:
  dy = collect("dy", path=path, info=False).squeeze()
  position = 0.5 * dy
  position[1:] += 0.5 * dy[:-1]
  position = np.cumsum(position)

  for var, ax in zip(varlist, axs):
      data = collect(var, path=path, tind=-1, info=False).squeeze()
      ax.plot(position, data, label=path)

for var, ax in zip(varlist, axs):
    ax.set_ylabel(var)

plt.legend()
plt.tight_layout()
plt.show()
