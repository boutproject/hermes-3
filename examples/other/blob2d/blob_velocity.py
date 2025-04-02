import numpy as np


def blob_velocity(n, method="COM", return_index=False):
    """
    Calculate blob velocity in normalized time and normalized grid spacing

    Input: Blob density as a 3D vector in the form  n[t,x,z] where t is time and x,z are the perpendicular spatial coordinates

    Keywords
    --------

    method    Method to use:
              'peak' -> Calculate velocity of the peak density
              'COM' -> Calculate centre of mass velocity
    return_index   return indices used to create velocity

    """
    from boututils import calculus as Calc

    size = n.shape

    x = np.zeros(size[0])
    z = np.zeros(size[0])

    if method == "peak":
        for i in np.arange(size[0]):
            nmax, nmin = np.amax((n[i, :, :])), np.amin((n[i, :, :]))
            xpos, zpos = np.where(n[i, :, :] == nmax)
            x[i] = xpos[0]
            z[i] = zpos[0]

    elif method == "COM":
        x = np.zeros(size[0])
        z = np.zeros(size[0])
        for i in np.arange(size[0]):
            data = n[i, :, :] - n[0, 0, 0]  # use corner cell rather than nmin
            ntot = np.sum(data[:, :])

            z[i] = np.sum(np.sum(data[:, :], axis=0) * (np.arange(size[2]))) / ntot
            x[i] = np.sum(np.sum(data[:, :], axis=1) * (np.arange(size[1]))) / ntot

    else:
        raise ValueError("Invalid method {}".format(method))

    vx = Calc.deriv(x)
    vz = Calc.deriv(z)

    if return_index:
        return vx, vz, x, z
    else:
        return vx, vz


def calculate_SI(path, method="COM"):
    from boutdata.collect import collect

    n = collect("Ne", path=path, info=False)

    vx, vy, xx, yy = blob_velocity(n[:, 2:-2, 0, :], method=method, return_index=True)

    # Calculate normalisation
    # Values are now in cells per output step
    dx = collect("dx", path=path, info=False)[0, 0]
    dz = collect("dz", path=path, info=False)[0, 0]

    g_11 = collect("g_11", path=path, info=False)[0, 0]
    g_33 = collect("g_33", path=path, info=False)[0, 0]

    wci = collect("Omega_ci", path=path, info=False)
    rho_s0 = collect("rho_s0", path=path, info=False)

    t_array = collect("t_array", path=path, info=False)
    dt = t_array[1] - t_array[0]

    dx_SI = dx * np.sqrt(g_11) * rho_s0  # X grid spacing in meters
    dz_SI = dz * np.sqrt(g_33) * rho_s0  # Z grid spacing in meters
    dt_SI = dt / wci  # Time between steps in seconds

    return vx * dx_SI / dt_SI, vy * dz_SI / dt_SI, xx, yy


if __name__ == "__main__":
    import pickle
    import sys

    path = sys.argv[1]  # First argument is the path to the data

    vx, vz, xx, yy = calculate_SI(path)

    f = open("Velocity.dat", "wb")
    pickle.dump(vx, f)
    f.close()

    f = open("Position.dat", "wb")
    pickle.dump(xx, f)
    f.close()

    try:
        import matplotlib.pyplot as plt

        plt.plot(vx)
        plt.show()
    except ImportError:
        pass
