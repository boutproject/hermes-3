import xhermes
from pathlib import Path
import matplotlib.pyplot as plt
import shutil
from pathlib import Path
import xhermes
import tempfile

# Look in the folder where this script is located for the simulation data
here = Path(__file__).parent

# Open the simulation data using xhermes
with tempfile.TemporaryDirectory() as tmpdirname:
    dst = Path(tmpdirname) / Path(here).name
    print(f"Copying {here} to tmpdir {dst}")
    shutil.copytree(src=here, dst=dst)
    ds = xhermes.open(dst)

fig, axs = plt.subplots(nrows=2)
(1.0 - ds["detachment_front_location"]).plot(ax=axs[0])
ds["detachment_control_src_mult"].plot(ax=axs[1])

axs[0].set_ylabel("error")
axs[1].set_ylabel("response")

plt.suptitle("Detachment control")
plt.tight_layout()

# # Save the figure into the folder where this script is located.
plt.savefig(here / "control.png")
