from boutdata import collect
import matplotlib.pyplot as plt
import numpy as np
import os

Nnorm = collect("Nnorm")
Tnorm = collect("Tnorm")

species = ["e", "d", "d+", "ne", "ne+", "ne+2", "ne+3", "ne+4"]
ncolor = ['b', 'b', 'b', (0,0,0.4), (0,0,0.5), (0,0,0.6), (0,0,0.8), (0,0,1)]
tcolor = ['r', 'r', 'r', (0.4,0,0), (0.5,0,0), (0.6,0,0), (0.8,0,0), (1,0,0)]
style = ['--', ':', '-.', '-', '-', '-', '-', '-']

skip = 10

N = [collect("N" + name)[::skip,0,:,0] for name in species]
T = [collect("P" + name)[::skip,0,:,0] / np.clip(n, 1e-5, None)
     for name, n in zip(species, N)]

n_max = max([np.amax(n) for n in N]) * Nnorm
t_max = max([np.amax(t) for t in T]) * Tnorm

for i in range(N[0].shape[0]):
    fig, ax1 = plt.subplots()
    ax1.set_xlabel(r'Cell number')
    ax2 = ax1.twinx()

    ax1, ax2 = ax2, ax1 # Swap left-right
    
    ax1.set_ylabel('Density [m^-3]', color='b')
    ax1.tick_params('y', colors='b')
    
    ax2.set_ylabel("Temperature [eV]", color='r')
    ax2.tick_params('y', colors='r')

    ax1.set_ylim(bottom=0.0, top=n_max)
    ax2.set_ylim(bottom=0.0, top=t_max)

    for name, n, c, s in zip(species, N, ncolor, style):
        ax1.plot(n[i] * Nnorm, label=name.capitalize(), color=c, linestyle=s)

    for name, t, c, s in zip(species, T, tcolor, style):
        ax2.plot(t[i] * Tnorm, label=name.capitalize(), color=c, linestyle=s)

    ax1.legend(loc="upper right")
    ax2.legend(loc="upper left")
    plt.show()
    plt.title(f"Time {i}")
    plt.savefig(f"1d_neon_{i:02d}.png")
    plt.close()

os.system('convert -resize 50% -delay 20 -loop 0 *.png 1d_neon.gif')
