Linear device simulations
=========================

Simulations without magnetic curvature. The diamagnetic current
term is usually turned off, because it has zero divergence.

Annulus simulations
-------------------

These set up a simulation where the Z coordinate is the azimuthal angle,
and X is a radial coordinate. This coordinate system has a pole
at the axis, and so X covers a range of radius with a hole in the middle.

annulus-isothermal-he
~~~~~~~~~~~~~~~~~~~~~

Pure isothermal helium plasma. No neutrals or atomic physics.

annulus-te-he
~~~~~~~~~~~~~

Pure helium plasma. No neutrals or atomic physics. The electron temperature
is evolved, but the ion temperature is fixed

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
