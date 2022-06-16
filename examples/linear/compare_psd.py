# Calculate spatial correlation


import numpy as np
import matplotlib.pyplot as plt
import pickle
import os

datasets = [("D 64x64", "annulus-isothermal-d", 'b--', 1, 1),
            ("D 128x128", "annulus-isothermal-d-2", 'b-', 1, 1),
            ("He", "annulus-isothermal-he", 'r--', 1, 0.5),
            (r"He $f\times\sqrt{2}$", "annulus-isothermal-he", 'r--', np.sqrt(2), 1)]


for label, path, style, fmult, alpha in datasets:
    with open(os.path.join(path, "psd.pkl"), "rb") as f:
        freq = pickle.load(f)
        psd = pickle.load(f)

    plt.plot(freq[1:] * fmult, psd[1:], style, label=label, alpha=alpha)

plt.yscale('log')
plt.xscale('log')
plt.xlabel('Frequency [kHz]')
plt.ylabel('Power Spectral Density')
plt.ylim([1e30, 1e36])
plt.legend()

# 17m device, half wavelength. 5eV.
# Isothermal sound wave speed
length = 17 # m
Te = 5 # eV
he_freq = 2 * np.sqrt(Te * 1.602e-19 / (4.*1.67e-27)) / length / 1000 # khz
d_freq = 2 * np.sqrt(Te * 1.602e-19 / (2.*1.67e-27)) / length / 1000

plt.savefig("psd-comparison-loglog.png")
plt.savefig("psd-comparison-loglog.pdf")

plt.xscale('linear')
plt.xlim([0,20])

plt.savefig("psd-comparison-loglin.png")
plt.savefig("psd-comparison-loglin.pdf")

plt.show()

