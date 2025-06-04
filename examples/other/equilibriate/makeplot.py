from boutdata import collect
import matplotlib.pyplot as plt

path = "."
names = ["Te", "Td+"]

# Read time axis
wci = collect("Omega_ci", path=path)
time = collect("t_array", path=path) * 1e3 / wci  # ms

# Read and plot temperatures
Tnorm = collect("Tnorm", path=path)
for name in names:
    T = collect(name, path=path).squeeze() * Tnorm
    plt.plot(time, T, label=name)

plt.legend()
plt.xlabel("Time [ms]")
plt.ylabel("Temperature [eV]")

plt.savefig("equilibriate.png")
plt.savefig("equilibriate.pdf")

plt.show()
