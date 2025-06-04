from boutdata import collect
import matplotlib.pyplot as plt
import numpy as np
import os

# Electron pressure
pi = collect("Pi")[::2,0,:,0]

yrange = [np.amin(pi), np.amax(pi)]

for i in range(pi.shape[0]):
    plt.plot(pi[i])
    plt.title(f"Time {i}")
    plt.ylim(*yrange)
    plt.xlabel("x index");
    plt.ylabel("Ion pressure");
    plt.savefig(f"1d_te_ti_{i:02d}.png")
    plt.close()

os.system('convert -resize 50% -delay 20 -loop 0 *.png 1d_te_ti.gif')
