#!/usr/bin/env python3
#
# Plots logarithm of time derivatives as a function of output step
#

from boutdata import collect
import numpy as np
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(description="Plot logarithm of time derivatives")
parser.add_argument("datapath", type=str, help="The path to the data (BOUT.dmp files)")

args = parser.parse_args()
path = args.datapath

varlist = ["NVd", "NVd+", "Nd", "Nd+", "Pd", "Pd+", "Pe"]

wci = float(collect("Omega_ci", path=path))

try:
    scale_timederivs = collect("scale_timederivs", path=path)[:,2:-2,:,0]
except ValueError:
    scale_timederivs = 1

for var in varlist:
    name = "ddt("+var+")"
    try:
        ddt = collect(name, path=path)[:,2:-2,:,0] / scale_timederivs
        rms = np.sqrt(np.mean(ddt**2, axis=(1,2))) * wci
        plt.plot(rms, label=name)
    except ValueError:
        print("Couldn't find " + name)

plt.legend()
plt.yscale('log')
plt.xlabel('Time steps')
plt.ylabel('Time derivative [Norm / s]')
plt.show()

