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

.. math::
   
   c_s = \sqrt{\sum_i P_i / \sum_i m_in_i}

This is set in the state as `sound_speed`, and is used for the numerical
diffusion terms in the parallel advection.

collisions
~~~~~~~~~~

For collisions between charged particles. In the following all quantities are
in SI units except the temperatures: `T` is in eV, so `eT` has units of Joules.

Debye length `\lambda_D`

.. math::

   \lambda_D = \sqrt{\frac{\epsilon_0 T_e}{n_e e}}
   
Coulomb logarithm, from [NRL formulary 2019], adapted to SI units

- For thermal electron-electron collisions

  .. math::

     \ln \lambda_ee = 16.6 - \frac{1}{2} \ln\left(n_e\right) + \frac{5}{4}\ln\left(T_e\right) - sqrt{10^{-5} + \left(\ln T_e - 2\right)^2 / 16} 

  
- Electron-ion collisions

  .. math::

     \ln \lambda_{ei} = \left\{\begin{array}{ll}
                              16.1 - \frac{1}{2}\ln\left(n_e\right) - \ln(Z) + \frac{3}{2}\ln\left(T_e\right) & \textrm{if} T_im_e/m_i < T_e < 10Z^2 \\
                              17.1 - \frac{1}{2}\ln\left(n_e\right) + \ln\left(T_e\right) & \textrm{if} T_im_e/m_i < 10Z^2 < T_e \\
                              9.09 - \frac{1}{2}\ln\left(n_i\right) + \frac{3}{2}\ln\left(T_i\right) - \ln\left(Z^2\mu\right) & \textrm{if} T_e < T_im_e/m_i \\
                              \end{array}\right.
     
- Mixed ion-ion collisions
  
  .. math::

     \ln \lambda_{ii'} = 16.1 - ln\left[\frac{ZZ'\left(\mu + \mu'\right)}{\mu T_{i'} + \mu'T_i}\left(\frac{n_iZ^2}{T_i} + \frac{n_{i'} Z'^2}{T_{i'}}\right)^{1/2}\right]


The frequency of charged species `a` colliding with charged species `b` is

.. math::

   \nu_{ab} = \frac{1}{3\pi^{3/2}\epsilon_0^2}\frac{Z_a^2 Z_b^2 n_b \ln\Lambda}{\left(v_a^2 + v_b^2\right)^{3/2}}\frac{\left(1 + m_a / m_b\right)}{m_a^2}


Note that the cgs expression in Hinton is divided by `\left(4\pi\epsilon_0\right)^2` to get
the expression in SI units.

For conservation of momentum, the collision frequencies `\nu_{ab}` and `\nu_{ba}` are
related by:

.. math::

   m_a n_a \nu_{ab} = m_b n_b \nu_{ba}

Momentum exchange, force on species `a` due to collisions with species `b`:

.. math::

   F_{ab} = \nu_{ab} m_a n_a \left( u_b - u_a \right)

   
Energy exchange, heat transferred to species `a` from species `b`:

.. math::

   Q_{ab} = \nu_{ab}\frac{3n_a m_a\left(T_b - T_a\right)}{m_a + m_b}

- Neutral-neutral collisions

  The cross-section is given by

  .. math::
     
     \sigma = \pi \left(\frac{d_1 + d_2}{2}\right)^2

  where `d_1` and `d_2` are the kinetic diameters of the two species. Typical values
  are [Wikipedia] for H2  2.89e-10m, He  2.60e-10m, Ne 2.75e-10m. 

  The mean relative velocity of the two species is

  .. math::

     v_{rel} = \sqrt{\frac{eT_1}{m_1} + \frac{eT_2}{m_2}}

  and so the collision rate of species 1 on species 2 is:

  .. math::

  \nu_{12} = v_{rel} n_2 \sigma


recycling
~~~~~~~~~

This component calculates the flux of a species into a Y boundary,
due to recycling of flow out of the boundary of another species.

The boundary fluxes might be set by sheath boundary conditions,
which potentially depend on the density and temperature of all species.
Recycling therefore can't be calculated until all species boundary conditions
have been set. It is therefore expected that this component is a top-level
component which comes after boundary conditions are set.


Atomic and molecular reactions
------------------------------

The formula for the reaction is used as the name of the component. This
makes writing the input file harder, since the formula must be in the exact same format
(e.g. `h + e` and `e + h` won't be recognised as being the same thing),
but makes reading and understanding the file easier.

To include a set of reactions, it is probably easiest to group them,
and then include the group name in the components list::

  [hermes]
  components = ..., reactions

  [reactions]
  type = (
          h + e -> h+ + 2e,  # ionisation
          h+ + e -> h,    # Radiative + 3-body recombination
         )

Note that brackets can be used to split the list of reactions over multiple lines,
and trailing commas are ignored. Comments can be used if needed to add explanation.
The name of the section does not need to be `reactions`, and multiple components could
be created with different reaction sets. Be careful not to include the same reaction
twice.

When reactions are added, all the species involved must be included, or an exception
should be thrown.

Hydrogenic processes
~~~~~~~~~~~~~~~~~~~~

Multiple isotopes of hydrogen can be evolved, so to keep track of this the
species labels `h`, `d` and `t` are all handled by the same hydrogen atomic
rates calculation. The following might therefore be used::

  [hermes]
  components = d, t, reactions

  [reactions]
  type = (
          d + e -> d+ + 2e,  # Deuterium ionisation
          t + e -> t+ + 2e,  # Tritium ionisation
         )

+------------------+-------------------------------------+
| Reaction         | Description                         |
+==================+=====================================+
| h + e -> h+ + 2e | Hydrogen ionisation (Amjuel 2.1.5)  |
| d + e -> d+ + 2e | Deuterium ionisation (Amjuel 2.1.5) |
| t + e -> t+ + 2e | Tritium ionisation (Amjuel 2.1.5)   |
+------------------+-------------------------------------+

Helium
~~~~~~

+----------------------+------------------------------------------------------------+
| Reaction             | Description                                                |
+======================+============================================================+
| he + e -> he+ + 2e   | He ionisation, unresolved metastables (Amjuel 2.3.9a)      |
| he+ + e -> he        | He+ recombination, unresolved metastables (Amjuel 2.3.13a) |
+----------------------+------------------------------------------------------------+

Neon
~~~~

These rates are taken from ADAS (96): SCD and PLT are used for the ionisation
rate and radiation energy loss; ACD and PRB for the recombination rate and radiation
energy loss. The ionisation potential is also included as a source or sink of energy
for the electrons.

+------------------------+-------------------------------------+
| Reaction               | Description                         |
+========================+=====================================+
| ne + e -> ne+ + 2e     | Neon ionisation                     |
| ne+ + e -> ne+2 + 2e   |                                     |
| ne+2 + e -> ne+3 + 2e  |                                     |
| ne+3 + e -> ne+4 + 2e  |                                     |
| ne+4 + e -> ne+5 + 2e  |                                     |
| ne+5 + e -> ne+6 + 2e  |                                     |
| ne+6 + e -> ne+7 + 2e  |                                     |
| ne+7 + e -> ne+8 + 2e  |                                     |
| ne+8 + e -> ne+9 + 2e  |                                     |
| ne+9 + e -> ne+10 + 2e |                                     |
| ne+ + e -> ne          | Neon recombination                  |
| ne+2 + e -> ne+        |                                     |
| ne+3 + e -> ne+2       |                                     |
| ne+4 + e -> ne+3       |                                     |
| ne+5 + e -> ne+4       |                                     |
| ne+6 + e -> ne+5       |                                     |
| ne+7 + e -> ne+6       |                                     |
| ne+8 + e -> ne+7       |                                     |
| ne+9 + e -> ne+8       |                                     |
| ne+10 + e -> ne+9      |                                     |
+------------------------+-------------------------------------+

