#!/usr/bin/env python3
#
# Plot a time-history of a solve, across different resolutions

import matplotlib.pyplot as plt
import numpy as np

from boutdata import collect

history = [
    # This is a list with one tuple per mesh used
    ("start",
     # Runs that led to the next mesh
     # Each is a directory containing the BOUT.dmp files
     ["."],
     # Runs not used for next mesh
     []),
]

varname = "Nd+"

grid_start = 0
last_grid_rms = None
for grid, paths, paths2 in history:

    wci = float(collect("Omega_ci", path=paths[0]))
    last_walltime = grid_start
    grid_first = True

    for path in paths:
        wtime = collect("wtime", path=path)
        wtime[0] += last_walltime
        wall_time = np.cumsum(wtime)
        last_walltime = wall_time[-1]
        try:
            scale_timederivs = collect("scale_timederivs", path=path)[:,2:-2,:,0]
            linestyle = '--'
        except ValueError:
            linestyle = '-'
            scale_timederivs = 1

        ddt = collect("ddt(" + varname + ")", path=path)[:,2:-2,:,0] / scale_timederivs
        rms = np.sqrt(np.mean(ddt**2, axis=(1,2))) * wci

        plt.plot(wall_time / (60 * 60), rms, linestyle=linestyle, color='k')
        if grid_first:
            grid_first = False
            plt.text(grid_start / (60*60) + 0.1, rms[0], grid)
            if last_grid_rms is not None:
                plt.plot([grid_start / (60 * 60)]*2, [last_grid_rms, rms[0]], 'r-o')
                plt.arrow(grid_start / (60 * 60), np.sqrt(last_grid_rms * rms[0]),
                          0, 0.1*np.sqrt(rms[0] * last_grid_rms),
                          length_includes_head=True, color='r',
                          head_width = 0.1, head_length = 0.2*np.sqrt(rms[0] * last_grid_rms))
    last_grid_rms = rms[-1]
    grid_start = last_walltime

    for path in paths2:
        wtime = collect("wtime", path=path)
        wtime[0] += last_walltime
        wall_time = np.cumsum(wtime)
        last_walltime = wall_time[-1]
        try:
            scale_timederivs = collect("scale_timederivs", path=path)[:,2:-2,:,0]
            linestyle = '--'
        except ValueError:
            linestyle = '-'
            scale_timederivs = 1

        ddt = collect("ddt(" + varname + ")", path=path)[:,2:-2,:,0] / scale_timederivs
        rms = np.sqrt(np.mean(ddt**2, axis=(1,2))) * wci

        plt.plot(wall_time / (60 * 60), rms, linestyle=linestyle, color='k')
    plt.plot([last_walltime / (60 * 60)], [rms[-1]], 'b*', markersize=10)

plt.yscale('log')
plt.xlabel("Wall time [hours]")
plt.ylabel(varname + r' time derivative [$10^{19}$m$^{-3}$ / sec]')
plt.savefig(f"history.png")
plt.savefig(f"history.pdf")
plt.show()
