from boutdata import collect
import matplotlib.pyplot as plt
import numpy as np
import os

Nnorm = collect("Nnorm")
Tnorm = collect("Tnorm")

#species = ["e", "d", "d+", "ne", "ne+2", "ne+4", "ne+6", "ne+8", "ne+10"]
#ncolor = ['b', 'b', 'b', (0,0,0.4), (0,0,0.5), (0,0,0.6), (0,0,0.7), (0,0,0.8), (0,0,1)]
#tcolor = ['r', 'r', 'r', (0.4,0,0), (0.5,0,0), (0.6,0,0), (0.7,0,0), (0.8,0,0), (1,0,0)]
#style = ['--', ':', '-.', '-', '-', '-', '-', '-', '-']

species = ["e", "d", "d+", "ne+6", "ne+8", "ne+10"]
ncolor = ['b', 'b', 'b', (0,0,0.4), (0,0,0.65), (0,0,1)]
tcolor = ['r', 'r', 'r', (0.4,0,0), (0.65,0,0), (1,0,0)]
style = ['--', ':', '-.', '-', '-', '-']

skip = 20

nfloor = 1e-5

dy = collect("dy")[0,:]
print(dy.shape)
ypos = np.ndarray((len(dy),))
ypos[0] = 0.5 * dy[0]
for i in range(1,len(dy)):
    ypos[i] = ypos[i-1] + 0.5*(dy[i-1] + dy[i])

N = [collect("N" + name)[::skip,0,:,0] for name in species]
T = [collect("P" + name)[::skip,0,:,0] / np.clip(n, nfloor, None)
     for name, n in zip(species, N)]

n_max = max([np.amax(n) for n in N]) * Nnorm
t_max = max([np.amax(t) for t in T]) * Tnorm

for i in range(N[0].shape[0]):
    fig, ax1 = plt.subplots(figsize=(10,5))
    ax1.set_xlabel(r'Parallel distance [m]')
    #ax1.set_xlim(left=-5, right=35)
    ax2 = ax1.twinx()

    ax1, ax2 = ax2, ax1 # Swap left-right
    
    ax1.set_ylabel('Density [m^-3]', color='b')
    ax1.tick_params('y', colors='b')
    ax1.set_yscale('log')

    ax2.set_ylabel("Temperature [eV]", color='r')
    ax2.tick_params('y', colors='r')

    ax1.set_ylim(bottom=nfloor * Nnorm, top=n_max)
    ax2.set_ylim(bottom=0.0, top=t_max)

    for name, n, c, s in zip(species, N, ncolor, style):
        ax1.plot(ypos, n[i] * Nnorm, label=name.capitalize(), color=c, linestyle=s)
        ax1.text(0.0, n[i][0] * Nnorm, name.capitalize(), color='b')

    for name, n, t, c, s in zip(species, N, T, tcolor, style):
        inds = np.argwhere(n[i] > nfloor) # Only show temperature when density sufficiently high
        ax2.plot(ypos[inds], t[i][inds] * Tnorm, label=name.capitalize(), color=c, linestyle=s)

    l1 = ax1.legend(bbox_to_anchor=(0.95, 0.0), loc="lower right", ncol=3)
    l2 = ax2.legend(bbox_to_anchor=(0.05, 0.0), loc="lower left", ncol=3)

    #ax1.set_title(f"Time {i}")

    fig.savefig(f"1d_neon_{i:02d}.png", bbox_extra_artists=(l1, l2))
    fig.savefig(f"1d_neon_{i:02d}.pdf", bbox_extra_artists=(l1, l2))

    plt.show()
    plt.close()

#os.system('convert -resize 50% -delay 20 -loop 0 *.png 1d_neon.gif')
