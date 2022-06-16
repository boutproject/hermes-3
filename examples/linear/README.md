Linear device simulations
=========================

Simulations without magnetic curvature. The diamagnetic current
term is usually turned off, because it has zero divergence.

Isothermal annulus simulations
------------------------------

These set up a simulation where the Z coordinate is the azimuthal angle,
and X is a radial coordinate. This coordinate system has a pole
at the axis, and so X covers a range of radius with a hole in the middle.

- annulus-isothermal-d: Pure isothermal deuterium plasma. No neutrals
  or atomic physics. Source of plasma ions uniform along the length of
  the device.

- annulus-isothermal-d-2: The same as `annulus-isothermal-d` but
  double resolution in directions perpendicular to the magnetic field.

- annulus-isothermal-he: Pure isothermal helium plasma. No neutrals or
  atomic physics.

Python scripts to analyse and compare simulations are stored in this
directory. The Power Spectral Density (PSD) from isothermal deuterium
and helium simulations is shown in the figure below:

![Power Spectral Density (PSD) comparison](psd-comparison-loglin.png)

Fixed neutrals annulus simulations
----------------------------------

These simulations provide a source of power to the electrons, which then
ionise a fixed background of neutrals to provide a plasma ion source.

- annulus-te-he-fixedneutrals: Evolves the electron temperature and
  helium ion (1+ only), with a fixed background of stationary helium
  atoms. Atomic reactions provide the density source for the plasma,
  so the only external input is the power.

Notes
-----

Useful routine for separating Z varying 'ac' components from 'dc' background:
```
def acdc(f):
    dc = np.mean(f, axis=-1)
    ac = np.copy(f)
    for z in range(ac.shape[-1]):
        ac[...,z] -= dc
    return ac, dc
```
