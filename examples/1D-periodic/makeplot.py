from boutdata import collect
import matplotlib.pyplot as plt

p = collect("Pi")

fig, ax = plt.subplots(1,1, figsize=(3,3))

for time in [0, 10, 25]:
    ax.plot(p[time, 0, :, 0])

plt.tight_layout()

plt.savefig("periodic1d.png")

plt.show()
