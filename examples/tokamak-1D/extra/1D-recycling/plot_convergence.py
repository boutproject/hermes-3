# Analyse convergence towards steady state

paths = ["cvode", "."]
labels = ["CVODE", "Backward-Euler"]
linestyles = ["-", "--"]

import matplotlib.pyplot as plt
import numpy as np

from boutdata import collect

cumulative_calls = []
times = []
source_mags = []

for path in paths:
    try:
        # Number of calls
        ncalls = collect("ncalls", path=path)  # number of calls each output

        cumulative_calls.append(np.cumsum(ncalls))   # Summed over outputs

        t_array = collect("t_array", path=path)
        wci = collect("Omega_ci", path=path)
        times.append(t_array / wci) # Seconds

        # density source
        m = collect("density_feedback_src_mult_d+", path=path)
        source_mags.append(np.abs(m + 1))
    except:
        raise

for calls, source, label, linestyle in zip(cumulative_calls, source_mags, labels, linestyles):
    plt.plot(calls, source, label=label, linestyle=linestyle)

plt.yscale('log')
plt.xlabel('RHS evaluations')
plt.ylabel('Source magnitude [arb]')
plt.legend()
plt.savefig("source_rhsevals.pdf")
plt.show()

for time, source, label, linestyle in zip(times, source_mags, labels, linestyles):
    plt.plot(time, source, label=label, linestyle=linestyle)
plt.yscale('log')
plt.xlabel('Time [s]')
plt.ylabel('Source magnitude [arb]')
plt.legend()
plt.savefig("source_time.pdf")
plt.show()

# Time derivatives

for path, calls, label, linestyle in zip(paths, cumulative_calls, labels, linestyles):
    for varname, color, varlabel in [(r"ddt(Nd+)", 'k', r'$\frac{\partial}{\partial t}n_{d+}$'),
                                     ("ddt(Pd+)", 'r', r'$\frac{\partial}{\partial t}p_{d+}$'),
                                     ("ddt(NVd+)",'b', r'$\frac{\partial}{\partial t}nv_{||d+}$')]:
        dv = collect(varname, path=path)
        wci = collect("Omega_ci", path=path)
        dv_rms = np.sqrt(np.sum(dv[:,0,:,0]**2, axis=1)) * wci
        plt.plot(calls, dv_rms, label=label + " " + varlabel, linestyle=linestyle, color=color)
plt.yscale('log')
plt.ylabel(r"RMS time derivative over domain [s$^{-1}$]")
plt.xlabel("RHS evaluations")
plt.legend()
plt.savefig("timederivs_rhsevals.pdf")
plt.show()

