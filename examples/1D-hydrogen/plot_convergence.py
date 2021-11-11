# Analyse convergence towards steady state

paths = [".", "qn"]
labels = ["NK", "QN"]

import matplotlib.pyplot as plt
import numpy as np

from boutdata import collect

cumulative_calls = []
times = []
source_mags = []

for path in paths:
    # Number of calls
    ncalls = collect("ncalls", path=path)  # number of calls each output

    cumulative_calls.append(np.cumsum(ncalls))   # Summed over outputs

    t_array = collect("t_array", path=path)
    wci = collect("Omega_ci", path=path)
    times.append(t_array / wci) # Seconds

    # density source
    m = collect("density_source_multiplier_d+", path=path)
    source_mags.append(np.abs(m + 1))

for calls, source, label in zip(cumulative_calls, source_mags, labels):
    plt.plot(calls, source, label=label)

plt.yscale('log')
plt.xlabel('RHS evaluations')
plt.ylabel('Source magnitude [arb]')
plt.legend()
plt.savefig("source_rhsevals.pdf")
plt.show()

for time, source, label in zip(times, source_mags, labels):
    plt.plot(time, source, label=label)
plt.yscale('log')
plt.xlabel('Time [s]')
plt.ylabel('Source magnitude [arb]')
plt.legend()
plt.savefig("source_time.pdf")
plt.show()

# Time derivatives

for path, calls, label in zip(paths, cumulative_calls, labels):
    for varname in ["ddt(Nd+)", "ddt(Pd+)", "ddt(NVd+)"]:
        dv = collect(varname, path=path)
        dv_rms = np.sqrt(np.sum(dv[:,0,:,0]**2, axis=1))
        plt.plot(calls, dv_rms, label=label + " " + varname)
plt.yscale('log')
plt.ylabel("RMS time derivative over domain")
plt.xlabel("RHS evaluations")
plt.legend()
plt.savefig("timederivs_rhsevals.pdf")
plt.show()

