# Electromagnetic Alfven wave with finite electron mass

![Alfven wave speed](alfven-wave.png)

## Equations

The ion density is fixed, with zero velocity and fixed temperature:
```
[i]
type = fixed_density, fixed_velocity, fixed_temperature

charge = 1
AA = 1

density = 1e19 # [m^-3]
velocity = 0
temperature = 100  # eV
```

The electron density is set to the ion density by quasi-neutrality,
and only the parallel momentum is evolved:
```
[e]
type = quasineutral, evolve_momentum, fixed_temperature

charge = -1
AA = 1./1836

temperature = 100 # eV
```
Note that the momentum is the canonical momentum because the
`electromagnetic` component calculates and corrects for the vector
potential term in the parallel momentum of all (charged) species.

Finally, the potential is evolved by a vorticity equation.

## Slab domain

The domain is a thin slab, with one cell in the X (radial) direction,
and periodic in both Y (parallel) and Z (binormal) directions.
The Y direction is used to set k_parallel, and the Z direction sets k_perp.

The number of cells in each dimension are specified by `nx`, `ny` and
`nz`. Note that `nx` includes 2 boundary cells on either side, so `nx
= 5` is one cell in the middle, and 4 boundary cells.

The size of the domain in each dimension in meters is set by `Lx`, `Ly` and `Lz`.
These are not directly used by Hermes-3 or BOUT++, but are used to calculate the
metric tensor components and cell sizes.

The Y and Z cell sizes are `dy = Ly / ny` and `dz = Lz / nz`.
The X cell size is more complicated because the field-aligned
operators assume a Clebsch coordinate system, in which B = nabla z
cross nabla x.  See the [coordinates manual page](https://bout-dev.readthedocs.io/en/stable/user_docs/coordinates.html#magnetic-field)
for details. The `dx` cell size therefore includes a factor of the
magnetic field strength `B`. With this choice the metric tensor is
diagonal, with elements
```
g11 = B^2
g22 = 1
g33 = 1
```
and the Jacobian is `J = 1 / B`

## Boundary conditions

The radial boundary conditions on potential are Neumann for the
non-constant (AC) components, and zero for the constant (DC) component:
```
[vorticity:laplacian]
inner_boundary_flags = 2
outer_boundary_flags = 2
```
All other boundary conditions are Neumann.
