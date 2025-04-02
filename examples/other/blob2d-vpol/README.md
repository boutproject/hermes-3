Blob2D-vpol
===========

A seeded plasma filament in 2D. This version is isothermal and cold ion,
so only the electron density and vorticity are evolved. A sheath-connected
closure is used for the parallel current.

Includes the (electron) polarisation drift.

![Electron density Ne at three times, showing propagation to the right](blob2d.png)

The model components are
```
[hermes]
components = e, vorticity, sheath_closure, polarisation_drift
```

The electron component consists of two types:
```
[e]  # Electrons
type = evolve_ne, isothermal
```

The `evolve_ne` component type evolves the electron density `Ne`. This component
has several options, which are set in the same section e.g.
```
poloidal_flows = false  # Y flows due to ExB
```

The `isothermal` component type sets the temperature to be a constant, and using
the density then sets the pressure. The constant temperature is also
set in this `[e]` section:
```
temperature = 5  # Temperature in eV
```

The `vorticity` component uses the pressure to calculate the diamagnetic current,
so must come after the `e` component. This component then calculates the potential.
Options to control the vorticity component are set in the `[vorticity]` section.

The `sheath_closure` component uses the potential, so must come after `vorticity`.
Options are also set as
```
[sheath_closure]
connection_length = 10 # meters
```

The `polarisation_drift` component approximates the polarisation drift as a potential flow, 
and applies it to all charged species, including electrons.

Analysis
--------

To plot the blob velocity at the three times shown above, run:
```
$ python3 makeplots.py
```
