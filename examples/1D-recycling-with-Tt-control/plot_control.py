import xhermes
from pathlib import Path
import matplotlib.pyplot as plt
from boutdata.data import BoutOptionsFile

# Look in the folder where this script is located for the simulation data
here = Path(__file__).parent

# Open the simulation data using xhermes
ds = xhermes.open(here)
# Select the first 100 time slices where the system is evolving rapidly
ds = ds.isel(t=slice(None, 100))
# Open the BOUT.inp file using boutdata
options = BoutOptionsFile(here / "BOUT.inp")

fig, axs = plt.subplots(nrows=2, ncols=3, sharex=True)

ds["Ne"].isel(y=0).plot(ax=axs[0][0], _labels=False)
# Indicate the density_upstream setpoint using a black dashed line
axs[0][0].axhline(options["d+"]["density_upstream"], color="k", linestyle="--")
axs[0][0].set_title("Upstream value")
axs[0][0].set_ylabel("Density")

ds["Ne"].isel(y=-1).plot(ax=axs[0][1], _labels=False)
axs[0][1].set_title("Target value")

ds["density_feedback_src_mult_d+"].plot(ax=axs[0][2], _labels=False)
axs[0][2].set_title("Source multiplier")

ds["Te"].isel(y=0).plot(ax=axs[1][0], _labels=False, label="$T_e$")
ds["Td+"].isel(y=0).plot(ax=axs[1][0], _labels=False, label="$T_{d+}$")
axs[1][0].set_ylabel("Temperature")
axs[1][0].legend()

ds["Te"].isel(y=-1).plot(ax=axs[1][1], _labels=False)
# Indicate the temperature_setpoint using a black dashed line
axs[1][1].axhline(options["e"]["temperature_setpoint"], color="k", linestyle="--")
ds["Td+"].isel(y=-1).plot(ax=axs[1][1], _labels=False)
ds["temperature_feedback_src_mult_e"].plot(ax=axs[1][2], _labels=False)

axs[1][1].set_xlabel("Time")

plt.tight_layout()

# Save the figure into the folder where this script is located.
plt.savefig(here / "control.png")
