from boutdata import collect
import matplotlib.pyplot as plt
import numpy as np
import os

Nnorm = collect("Nnorm")
Tnorm = collect("Tnorm")

ni = collect("Nd+")[::10,0,:,0]  # Ion density
nn = collect("Nd")[::10,0,:,0]   # Neutral density
pe = collect("Pe")[::10,0,:,0]   # Electron pressure
pi = collect("Pd+")[::10,0,:,0]   # Ion pressure
pn = collect("Pd")[::10,0,:,0]   # Neutral pressure

Te = pe / np.clip(ni, 1e-5, None)  # Electron temperature
Ti = pi / np.clip(ni, 1e-5, None)  # Ion temperature
Tn = pn / np.clip(nn, 1e-5, None)  # Neutral temperature

n_max = max([np.amax(ni),np.amax(nn)]) * Nnorm
t_max = max([np.amax(Te),np.amax(Ti), np.amax(Tn)]) * Tnorm

for i in range(ni.shape[0]):
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

    ax1.plot(ni[i] * Nnorm, label="Ion density", color='b', linestyle='-')
    ax1.plot(nn[i] * Nnorm, label="Neutral density", color='b', linestyle='--')
    #ax1.legend("upper left")

    ax2.plot(Te[i] * Tnorm, label="Electron", color='r', linestyle=":")
    ax2.plot(Ti[i] * Tnorm, label="Ion", color='r', linestyle='-')
    ax2.plot(Tn[i] * Tnorm, label="Neutral", color='r', linestyle='--')
    ax2.legend(loc="upper right")
    
    plt.title(f"Time {i}")
    plt.savefig(f"1d_recycling_{i:02d}.png")
    plt.savefig(f"1d_recycling_{i:02d}.pdf")
    plt.close()

os.system('convert -resize 50% -delay 20 -loop 0 *.png 1d_recycling.gif')
