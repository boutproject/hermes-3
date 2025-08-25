Diffusion
=========

To run this simulation a grid file is needed. Assumed to be called `tokamak.nc`.


Solves simple diffusion perpendicular to the magnetic field,
of a density for a single species

    dn/dt = Div( D * Grad_perp(n))

This is done by specifying a single species:

```
[hermes]
components = h
```

and then specifying that the density should be evolved, and there
should be anomalous diffusion:
```
[h]
type = evolve_density, anomalous_diffusion
```

The coefficient of diffusion is taken from the setting
```
[h]

anomalous_D = 2    # Density diffusion [m^2/s]
```

The evolving density is `Nh`, created by prepending `N` to the species name.
The initial condition is flat (1 everywhere), and the
boundary conditions are set to high density on the core and low density on the edge
```
[Nh]

function = 1 

bndry_core = dirichlet(1.0)
bndry_pf = dirichlet(0.1)
bndry_sol = dirichlet(0.1)
```
