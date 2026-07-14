1D periodic domain, Te and Ti
=============================

A fluid is evolved in 1D, imposing quasineutrality and zero net current.
Both electron and ion pressures are evolved, but there is no exchange
of energy between them, or heat conduction.

The model components are ions (i), electrons (e), and a constraint
that the net current is zero. This constraint applies the electron
pressure to the ion momentum equation, and sets the electron velocity
to be equal to the ion velocity.
```
[hermes]
components = i, e, zero_current
```

The ion density, pressure and momentum equations are evolved:
```
[i]  # Ions
type = evolve_density, evolve_pressure, evolve_momentum
```

The electron density is set to the ion density by quasineutrality,
and only the electron pressure is evolved.
```
[e] # Electrons
type = quasineutral, evolve_pressure
```
