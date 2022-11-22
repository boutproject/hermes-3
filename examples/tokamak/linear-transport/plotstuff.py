import xhermes
import matplotlib.pyplot as plt

# Open the dataset
bd = xhermes.open(".", geometry="toroidal", gridfilepath="../tokamak.nc")

# Make a plot of electron temperature against time
bd["Te"].isel(x=30, theta=30, zeta=0).plot()
plt.yscale('log')
plt.tight_layout()
plt.savefig("Te_30_30.png")
plt.closeall()

bd["Te"].isel(t=100, zeta=0).bout.contourf()
plt.savefig("Te2D.png")
plt.closeall()

