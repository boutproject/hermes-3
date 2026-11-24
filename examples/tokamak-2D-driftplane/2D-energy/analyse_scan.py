#!/usr/bin/env python3
"""Plot total-energy histories and grid differences using xhermes inventories."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xhermes

cli = argparse.ArgumentParser(description=__doc__)
cli.add_argument("directory", type=Path)
args = cli.parse_args()
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
rows = []
for path in sorted(args.directory.glob("*x*"), key=lambda p: int(p.name.split("x")[0])):
    with xhermes.open(path, keep_xboundaries=True, keep_yboundaries=False) as ds:
        inventory = ds.hermes.inventory().compute()
        energy = inventory.ion_energy + inventory.electron_thermal
        change = (energy - energy.isel(t=0)) / energy.isel(t=0)
        axes[0].plot(ds.t, change, label=path.name)
        row = {
            "grid": path.name,
            "nx": int(path.name.split("x")[0]),
            "final_energy_J": float(energy.isel(t=-1)),
            "final_electron_thermal_J": float(inventory.electron_thermal.isel(t=-1)),
            "max_relative_drift": float(abs(change).max()),
        }
        row.update(
            {
                name: float(value.isel(t=-1))
                for name, value in inventory.items()
                if name.startswith("particles_")
            }
        )
        rows.append(row)
if len(rows) < 2:
    raise ValueError("At least two completed grids are needed for convergence")
frame = pd.DataFrame(rows)
difference = abs(np.diff(frame.final_electron_thermal_J))
axes[1].semilogx(frame.nx.iloc[:-1], difference, "o-")
axes[0].set(xlabel="Time [s]", ylabel="Relative energy change")
axes[0].legend()
axes[1].set(xlabel="Radial cells", ylabel="Adjacent electron-energy difference [J]")
fig.tight_layout()
fig.savefig(args.directory / "convergence.pdf")
frame.to_csv(args.directory / "convergence.csv", index=False)
print(frame.to_string(index=False))
