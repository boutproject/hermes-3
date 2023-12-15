import xhermes
import numpy as np


def extract_data(path):
    bd = xhermes.open(path)

    qe = 1.602e-19  # Coulombs
    eps0 = 8.854e-12  # F/m
    Mp = 1.67e-27  # kg

    # Electron temperature in eV
    Te = float(bd["Te"][2, 0, 0])
    Ti = float(bd["Ti"][2, 0, 0])
    # Electron density in m^-3
    Ne = bd.metadata["Nnorm"]
    # Ion charge
    Zi = bd.options["i"].evaluate_scalar("charge")
    Ni = Ne / Zi
    # Ion mass
    Mi = bd.options["i"].evaluate_scalar("AA") * Mp
    # Electron mass
    Me = bd.options["e"].evaluate_scalar("AA") * Mp

    # Perpendicular wavelength
    Lz = bd.options["mesh"].evaluate_scalar("Lz")  # meters
    kz = 2 * np.pi / Lz

    # Parallel wavelength
    Ly = bd.options["mesh"].evaluate_scalar("Ly")  # meters
    ky = 2 * np.pi / Ly

    # Inverse density scale length
    inv_Ln = bd.options["mesh"].evaluate_scalar("inv_Ln")  # 1/m

    # Magnetic field
    B = bd.options["mesh"].evaluate_scalar("B")

    # Coulomb Logarithm (electron-ion)
    CoulombLog = 31 - 0.5 * np.log(Ne) + np.log(Te)

    # Electron-ion collision time
    nu_ei = (
        qe**4
        * Zi**2
        * Ni
        * CoulombLog
        * (1 + Me / Mi)
        / (
            3
            * np.pi ** (3 / 2)
            * eps0**2
            * (2 * qe * (Te / Me + Ti / Mi)) ** (3 / 2)
            * Me**2
        )
    )

    wci = qe * B / Mi
    wce = qe * B / Me

    wstar = kz * Te * inv_Ln / B
    sigmapar = (ky / kz) ** 2 * wci * wce / (0.51 * nu_ei)
    sparn = sigmapar / wstar

    print(f"w* = {wstar} sigma_|| = {sigmapar}")
    print(f"sigma/w* = {sparn}")

    # Calculate analytic result
    # Note: This is for the simple resistive drift wave, as in Umansky 2008.
    #       It assumes zero electron mass.
    t = 0.5 * (np.sqrt(sparn**4 + 16 * sparn**2) - sparn**2)
    wr = 0.5 * np.sqrt(t)
    wi = sparn / np.sqrt(t) - 0.5 * sparn
    print(f"Simple analytic gamma / w* = {wi}, omega / w* = {wr}")

    # Cubic equation due to finite electron mass
    roots = np.roots([wstar / (0.51 * nu_ei), 1j, -sigmapar / wstar, sigmapar / wstar])
    # Find the root with the fastest growth rate
    w = roots[np.argmax(np.imag(roots))]
    wr = w.real  # Real frequency
    wi = w.imag  # Growth rate
    print(f"Analytic gamma / w* = {wi}, omega / w* = {wr}")

    # Calculate simulation growth rate and frequencies

    n = bd["Ni"][:, 2, 0, :].squeeze() - Ne  # Perturbation n(t,z)
    nrms = np.sqrt(np.mean(n**2, axis=-1))  # nrms(t)

    dt = float(bd["t"][1] - bd["t"][0])  # Seconds
    gamma = np.mean(np.gradient(np.log(nrms))[-10:]) / dt  # Growth rate [1/s]

    print(f"gamma = {gamma}, gamma / w* = {gamma / wstar}")

    # Track the motion of the peak to infer phase velocity
    n = np.array(n)
    nt, nz = n.shape
    peak = np.zeros(nt)

    for t in range(nt):
        ind = np.argmax(n[t, :])  # Index of maximum value

        # Refine maximum location by fitting a quadratic through 3 points
        c = n[t, ind]
        m = n[t, (ind - 1) % nz]
        p = n[t, (ind + 1) % nz]
        # Fit y = c + ax + bx**2
        a = 0.5 * (p - m)
        b = p - (c + a)
        peak[t] = ind - 0.5 * a / b  # Maximum

    # Check for periodic recurrence
    if peak[-1] > peak[-2]:
        # Increasing; watch for wrapping around the top
        for i in range(1, nt):
            if peak[i] < peak[i - 1]:
                peak[i:] += nz
    else:
        # Decreasing; watch for wrapping around the bottom
        for i in range(1, nt):
            if peak[i] > peak[i - 1]:
                peak[i:] -= nz

    omega = 2 * np.pi * np.mean(np.gradient(peak)[-10:]) / dt / nz
    print(f"omega = {omega}, w / w* = {omega / wstar}")

    return {
        "wstar": wstar,
        "sigmapar": sigmapar,
        "CoulombLog": CoulombLog,
        "nu_ei": nu_ei,
        "wci": wci,
        "wce": wce,
        "analytic_gamma": wi,
        "analytic_omega": wr,
        "gamma": gamma / wstar,
        "omega": omega / wstar,
    }
