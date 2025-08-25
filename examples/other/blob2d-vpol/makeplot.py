from boutdata import collect
import matplotlib.pyplot as plt

n = collect("Ne")

fig, axs = plt.subplots(1,3, figsize=(9,3))

for ax, time in zip(axs, [15, 30, 45]):
    ax.contourf(n[time, 2:-2, 0, :].T, 50)

plt.tight_layout()

plt.savefig("blob2d.png")

plt.show()
