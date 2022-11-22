import xhermes
import matplotlib.pyplot as plt
import numpy as np

# Open the dataset
bd = xhermes.open(".", geometry="toroidal", gridfilepath="eqd1031204007.00740.bout.nc")
Te = bd["Te"]

Te.isel(t=100, zeta=0).bout.pcolormesh()
plt.savefig("Te2D.png")
plt.close('all')

# Calculate decay rate (positive for decaying in time)
decay = np.log(Te[-2,:,:,0] / Te[-1,:,:,0]) / (Te.t[-1] - Te.t[-2])

for ix, iy in [(30, 30), (20, 50)]:
    # Make a plot of electron temperature against time
    Te.isel(x=ix, theta=iy, zeta=0).plot()

    gamma = float(decay[ix,iy])
    T1 = Te[-1,ix,iy,0]
    T0 = T1 / np.exp(-gamma * Te.t[-1])

    plt.plot([0.0, Te.t[-1]], [T0, T1], '--k')
    plt.title(r'Decay time: {:.2f} ms @ ix={},iy={}'.format(1000./gamma, ix, iy))
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig("Te_{}_{}.png".format(ix, iy))
    plt.close('all')
