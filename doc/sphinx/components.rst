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

Evolves the pressure in time. This pressure is named `P<species>` where `<species>`
is the short name of the evolving species e.g. `Pe`.

Species parallel dynamics
-------------------------

evolve_momentum
~~~~~~~~~~~~~~~

Evolves the momentum `NV<species>` in time. 

zero_current
~~~~~~~~~~~~

This imposes a zero-current Ohm's law, calculating a parallel
electric field which balances the electron pressure gradient.
This electric field is then used to calculate a force on the other species.

Collective quantities
---------------------

These components combine multiple species together

sound_speed
~~~~~~~~~~~

Calculates the collective sound speed, by summing the pressure of all species,
and dividing by the sum of the mass density of all species:

\[
c_s = \sqrt{\sum_i P_i / \sum_i m_in_i}
\]

This is set in the state as `sound_speed`, and is used for the numerical
diffusion terms in the parallel advection.

