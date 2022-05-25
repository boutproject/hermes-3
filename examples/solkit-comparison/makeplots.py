
from boutdata import collect

import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1,2, figsize=[12.8, 4.8])

colors = ['b', 'r', 'g', 'm']
col_iter = iter(colors)

# Single temperature case

for path in ["single-temperature"]:
  tnorm = collect("Tnorm", path=path)
  nnorm = collect("Nnorm", path=path)
  n_1 = collect("Nd+", path=path, tind=-1)[-1,0,:,0]
  pe_1 = collect("Pe", path=path, tind=-1)[-1,0,:,0]
  t_1 = pe_1 / n_1  # Plasma temperature
  c = next(col_iter)
  ax1.plot(t_1 * tnorm, color=c, label=path)
  ax2.plot(n_1 * nnorm, color=c, label=path)

# Two temperatures cases

for path in ["two-temperatures"]:
    n = collect("Nd+", path=path, tind=-1)[-1,0,:,0]
    pe = collect("Pe", path=path, tind=-1)[-1,0,:,0]
    pi = collect("Pd+", path=path, tind=-1)[-1,0,:,0]
    te = pe / n # Electron temperature
    ti = pi / n # Ion temperature

    c = next(col_iter)
    ax1.plot(te * tnorm, color=c, linestyle="--", label=r"{} $T_e$".format(path))
    ax1.plot(ti * tnorm, color=c, linestyle=":", label=r"{} $T_i$".format(path))
    ax1.plot(0.5 * (ti + te) * tnorm, color=c, linestyle="-", label=r"{} $(T_i + T_e)/2$".format(path))
    ax2.plot(n * nnorm, color=c, label=path)

ax1.set_ylabel("Temperature [eV]")
ax1.legend()

ax2.set_ylabel(r"Density [m$^{-3}$]")
ax2.legend()

plt.savefig("comparison.png")
plt.savefig("comparison.pdf")
plt.show()
