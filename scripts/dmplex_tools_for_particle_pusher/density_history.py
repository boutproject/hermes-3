#!/usr/bin/env python3

import xbout
import argparse
from pathlib import Path
import matplotlib.pyplot as plt


def plot(case_path):

    ds = xbout.load.open_boutdataset(
        datapath=case_path / "BOUT.dmp.*.nc",
        inputfilepath=case_path / "BOUT.inp",
        keep_xboundaries=False,
        keep_yboundaries=False,
        info=False,
    )

    fig, axes = plt.subplots(1, 3, figsize=(9, 3))

    dv = ds["dx"] * ds["dy"] * ds["dz"] * ds["J"]
    Nn_tot = (ds["neutral_density"] * dv).sum(("x", "y"))
    Nn_avg = ds["neutral_density"].mean(("x", "y"))
    Ni_avg = ds["ion_density"].mean(("x", "y"))

    ax = axes[0]
    Ni_avg.plot(ax=ax, marker="o")
    ax.set_ylabel("Normalised density")
    ax.set_title("Mean ion density")
    
    ax = axes[1]
    Nn_avg.plot(ax=ax, marker="o")
    ax.set_ylabel("Normalised density")
    ax.set_title("Mean particle density")

    ax = axes[2]
    ds["total_ion_mass"].plot(ax=ax, marker="o", c="coral", label="total_ion_mass")
    ds["total_neutral_mass"].plot(
        ax=ax, marker="d", c="skyblue", label="total_neutral_mass"
    )
    ds["total_mass"].plot(
        ax=ax, marker="v", c="black", alpha=0.3, lw=3, label="total_mass"
    )
    ax.legend(fontsize="x-small")
    ax.set_ylabel("-")
    ax.set_title("Density volume integrals")


    for ax in axes:
        ax.set_xlabel("time")

    fig.tight_layout()
    fig.savefig(f"density_history_{case_path.name}.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Density history plotter")
    parser.add_argument("case_path", type=str, help="Simulation directory")

    args = parser.parse_args()

    plot(Path(args.case_path))
