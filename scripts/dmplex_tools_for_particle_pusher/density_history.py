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

    fig, axes = plt.subplots(1,2, figsize = (6,3))

    dv = ds["dx"] * ds["dy"] * ds["dz"] * ds["J"]
    Ntot = (ds["neutral_density"] * dv).sum(("x", "y"))
    Navg = ds["neutral_density"].mean(("x", "y"))

    ax = axes[0]
    Ntot.plot(ax = ax, marker = "o")
    ax.set_xlabel("time")
    ax.set_ylabel("Total particle count")
    
    ax = axes[1]
    Navg.plot(ax = ax, marker = "o")
    ax.set_ylabel("Mean density")
    ax.set_title("Mean particle density")

    for ax in axes:
        ax.set_title("Particle count")

    fig.tight_layout()
    fig.savefig(f"density_history_{case_path.name}.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Density history plotter")
    parser.add_argument("case_path", type=str, help="Simulation directory")

    args = parser.parse_args()

    plot(Path(args.case_path))
