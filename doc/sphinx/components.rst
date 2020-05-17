.. _sec-components:

Components
==========



Species density
---------------

The density of a species can be calculated in several different ways,
and are usually needed by other components.

evolve_density
~~~~~~~~~~~~~~

This component evolves the species density in time, using the BOUT++
time integration solver.

fixed_fraction_ions
~~~~~~~~~~~~~~~~~~~

This sets the density of a species to a fraction of the electron density.

quasineutral
~~~~~~~~~~~~

This component sets the density of one species, so that the overall
charge density is zero everywhere. This must therefore be done after
all other charged species densities have been calculated. It only
makes sense to use this component for species with a non-zero charge.



Species pressure and temperature
--------------------------------

evolve_pressure
~~~~~~~~~~~~~~~

