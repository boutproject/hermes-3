from boutdata import collect
import matplotlib.pyplot as plt

path = "data"

n = collect("Ni", path=path)
nv = collect("NVi", path=path)
p = collect("Pi", path=path)

v = nv / n

time = collect("t_array", path=path)
wci = collect("Omega_ci", path=path)
time /= wci

nt = time.shape[0]

fig, (axn, axv, axp) = plt.subplots(1, 3)

for t in [0, nt // 2, nt - 1]:
    axn.plot(n[t,0,:,0], label=f"t={time[t]}s")
    axv.plot(v[t,0,:,0], label=f"t={time[t]}s")
    axp.plot(p[t,0,:,0], label=f"t={time[t]}s")

axn.set_title("Density")
axv.set_title("Velocity")
axp.set_title("Pressure")

plt.legend()
plt.savefig("toro-2.png")
    
plt.show()
