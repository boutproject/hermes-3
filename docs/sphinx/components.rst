.. _sec-components:

Components
==========

This section describes the model components currently available. 

Species density
---------------

The density of a species can be calculated in several different ways,
and are usually needed by other components.

.. _fixed_density:

fixed_density
~~~~~~~~~~~~~

Set the density to a value which does not change in time. For example:

.. code-block:: ini

   [d]
   type = fixed_density, ...

   AA = 2 # Atomic mass
   charge = 0
   density = 1e17 # In m^-3

Note that the density can be a function of `x`, `y` and `z` coordinates
for spatial variation.

The implementation is in the `FixedDensity` class:

.. doxygenstruct:: FixedDensity
   :members:

.. _evolve_density:

evolve_density
~~~~~~~~~~~~~~

This component evolves the species density in time, using the BOUT++
time integration solver.

The implementation is in the `EvolveDensity` class:

.. doxygenstruct:: EvolveDensity
   :members:

.. _upstream_density_feedback:

upstream_density_feedback
~~~~~~~~~~~~~~~~~~~~~~~~~

This is intended for 1D simulations, where the density at :math:`y=0` is set
by adjusting an input source. This component uses a PI controller method
to scale the density source up and down, to maintain the specified upstream
density.

For example:

.. code-block:: ini

   [d+]
   type = ..., upstream_density_feedback

   density_upstream = 1e19      # Density in m^-3
   density_controller_p = 1e-2  # Feedback controller proportional (p) parameter
   density_controller_i = 1e-3  # Feedback controller integral (i) parameter

   [Nd+]
   source = h(pi - y)  # Source shape

The implementation is in the `UpstreamDensityFeedback` class:

.. doxygenstruct:: UpstreamDensityFeedback
   :members:

fixed_fraction_ions
~~~~~~~~~~~~~~~~~~~

This sets the density of a species to a fraction of the electron density.

.. _quasineutral:

quasineutral
~~~~~~~~~~~~

This component sets the density of one species, so that the overall
charge density is zero everywhere. This must therefore be done after
all other charged species densities have been calculated. It only
makes sense to use this component for species with a non-zero charge.


Species pressure and temperature
--------------------------------

.. _isothermal:

isothermal
~~~~~~~~~~

Sets the temperature of a species to a fixed value

.. code-block:: ini

   [e]
   type = ..., isothermal

   temperature = 10   # Constant temperature [eV]

.. doxygenstruct:: Isothermal
   :members:

.. _evolve_pressure:

evolve_pressure
~~~~~~~~~~~~~~~

Evolves the pressure in time. This pressure is named `P<species>` where `<species>`
is the short name of the evolving species e.g. `Pe`.

By default parallel thermal conduction is included, which requires a collision
time. If collisions are not calculated, then thermal conduction should be turned off
by setting `thermal_conduction = false` in the input options.

Notes:

- Heat conduction through the boundary is turned off currently. This is because
  heat losses are usually calculated at the sheath, so any additional heat conduction
  would be in addition to the sheath heat transmission already included.

The implementation is in `EvolvePressure`:

.. doxygenstruct:: EvolvePressure
   :members:


Species parallel dynamics
-------------------------

.. _evolve_momentum:

evolve_momentum
~~~~~~~~~~~~~~~

Evolves the momentum `NV<species>` in time. The evolving quantity includes the atomic
mass number, so should be divided by `AA` to obtain the particle flux.

.. _zero_current:

zero_current
~~~~~~~~~~~~

This calculates the parallel flow of one charged species so that there is no net current,
using flows already calculated for other species. It is used like `quasineutral`:

.. code-block:: ini

   [hermes]
   components = h+, ..., e, ...   # Note: e after all other species
   
   [e]
   type = ..., zero_current,... # Set e:velocity

   charge = -1 # Species must have a charge


electron_force_balance
~~~~~~~~~~~~~~~~~~~~~~

This calculates a parallel electric field which balances the electron
pressure gradient and other forces on the electrons (including
collisional friction, thermal forces):

.. math::

   E_{||} = \left(-\nabla p_e + F\right) / n_e

where :math:`F` is the `momentum_source` for the electrons.
This electric field is then used to calculate a force on the other species:

.. math::

   F_z = Z n_z E_{||}

which is added to the ion's `momentum_source`. 

The implementation is in `ElectronForceBalance`:

.. doxygenstruct:: ElectronForceBalance
   :members:

Neutral gas models
------------------

The `neutral_mixed` component solves fluid equations along :math:`y`
(parallel to the magnetic field), and uses diffusive transport in :math:`x`
and :math:`z`.  It was adopted from the approach used in UEDGE and this paper
[Journal of Nuclear Materials, vol. 313-316, pp. 559-563 (2003)].

.. math::
   
   \begin{aligned}\frac{\partial n_n}{\partial t} =& -\nabla\cdot\left(n_n\mathbf{b}v_{||n} + n_n\mathbf{v}_{\perp n}\right) + S\\ \frac{\partial}{\partial t}\left(n_nv_{||n}\right) =& -\nabla\cdot\left(n_nv_{||n} \mathbf{b}v_{||n} + n_nv_{||n}\mathbf{v}_{\perp n}\right) - \partial_{||}p_n + \nabla_{||}\left(D_{nn}n_n\partial_{||}v_{||n}\right) + F \\ \frac{\partial p_n}{\partial t} =& -\nabla\cdot\left(p_n\mathbf{b}v_{||n} + p_n\mathbf{v}_{\perp n}\right) - \frac{2}{3}p_n\nabla\cdot\left(\mathbf{b}v_{||n}\right) + \nabla\cdot\left(D_{nn}n_n\nabla_\perp T_n\right) + \frac{2}{3}Q \end{aligned}

The parallel momentum is evolved, so that it can be exchanged with the
plasma parallel momentum, but the mass is neglected for perpendicular
motion. In the perpendicular direction, therefore, the motion is a
balance between the friction (primarily with the plasma through charge
exchange) and the pressure gradient:

.. math::

   \mathbf{v}_{\perp n} = -D_{nn}\frac{1}{p_n}\nabla_\perp p_n

At the moment there is no attempt to limit these velocities, which has
been found necessary in UEDGE to get physical results in better
agreement with kinetic neutral models [Discussion, T.Rognlien].

.. _noflow_boundary:

noflow_boundary
~~~~~~~~~~~~~~~

This is a species component which imposes a no-flow boundary condition
on y (parallel) boundaries.

- Zero-gradient boundary conditions are applied to `density`,
  `temperature` and `pressure` fields, if they are set.
- Zero-value boundary conditions are applied to `velocity` and
  `momentum` if they are set.

By default both yup and ydown boundaries are set, but can be turned
off by setting `noflow_lower_y` or `noflow_upper_y` to `false`.

Example: To set no-flow boundary condition on an ion `d+` at the lower
y boundary, with a sheath boundary at the upper y boundary:

.. code-block:: ini

   [hermes]
   components = d+, sheath_boundary

   [d+]
   type = noflow_boundary

   noflow_lower_y = true   # This is the default
   noflow_upper_y = false  # Turn off no-flow at upper y for d+ species

   [sheath_boundary]
   lower_y = false         # Turn off sheath lower boundary for all species
   upper_y = true

Note that currently `noflow_boundary` is set per-species, whereas
`sheath_boundary` is applied to all species. This is because sheath
boundary conditions couple all charged species together, and doesn't
affect neutral species.

The implementation is in `NoFlowBoundary`:

.. doxygenstruct:: NoFlowBoundary
   :members:

Collective quantities
---------------------

These components combine multiple species together. They are typically
listed after all the species groups in the component list, so that all
the species are present in the state.

One of the most important is the `collisions`_ component. This sets collision
times for all species, which are then used 

.. _sound_speed:

sound_speed
~~~~~~~~~~~

Calculates the collective sound speed, by summing the pressure of all species,
and dividing by the sum of the mass density of all species:

.. math::
   
   c_s = \sqrt{\sum_i P_i / \sum_i m_in_i}

This is set in the state as `sound_speed`, and is used for the numerical
diffusion terms in the parallel advection.

.. _neutral_parallel_diffusion:

neutral_parallel_diffusion
~~~~~~~~~~~~~~~~~~~~~~~~~~

This adds diffusion to **all** neutral species (those with no or zero charge),
because it needs to be calculated after the collision frequencies are known.
It is intended mainly for 1D simulations, to provide effective parallel
diffusion of particles, momentum and energy due to the projection of
cross-field diffusion:

.. math::

   \begin{aligned}
   \frac{\partial n_n}{\partial t} =& \ldots + \nabla\cdot\left(\mathbf{b}D_n n_n\partial_{||}p_n\right) \\
   \frac{\partial p_n}{\partial t} =& \ldots + \nabla\cdot\left(\mathbf{b}D_n p_n\partial_{||}p_n\right) + \frac{2}{3}\nabla\cdot\left(\mathbf{b}\kappa_n \partial_{||}T_n\right) \\
   \frac{\partial}{\partial t}\left(n_nv_{||n}\right) =& \ldots + \nabla\cdot\left(\mathbf{b}D_n n_nv_{||n} \partial_{||}p_n\right) + \nabla\cdot\left(\mathbf{b}\eta_n \partial_{||}T_n\right)
   \end{aligned}

The diffusion coefficient is calculated as

.. math::

   D_n = \left(\frac{B}{B_{pol}}\right)^2 \frac{T_n}{A \nu}

where `A` is the atomic mass number; :math:`\nu` is the collision
frequency. The factor :math:`B / B_{pol}` is the projection of the cross-field
direction on the parallel transport, and is the `dneut` input setting.

.. doxygenstruct:: NeutralParallelDiffusion
   :members:


.. _collisions:

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

where :math:`d_1` and :math:`d_2` are the kinetic diameters of the two
species. Typical values are [Wikipedia] for H2 2.89e-10m, He
2.60e-10m, Ne 2.75e-10m.

The mean relative velocity of the two species is

.. math::

   v_{rel} = \sqrt{\frac{eT_1}{m_1} + \frac{eT_2}{m_2}}

and so the collision rate of species 1 on species 2 is:

.. math::

   \nu_{12} = v_{rel} n_2 \sigma

The implementation is in `Collisions`:

.. doxygenstruct:: Collisions
   :members:

.. _thermal_force:

thermal_force
~~~~~~~~~~~~~

This implements simple expressions for the thermal force. If the
`electron_ion` option is true (which is the default), then a momentum
source is added to all ions:

.. math::

   F_z = 0.71 n_z Z^2 \nabla_{||}T_e

where :math:`n_z` is the density of the ions of charge :math:`Z`. There
is an equal and opposite force on the electrons.

If the `ion_ion` option is true (the default), then forces are
calculated between light species (atomic mass < 4) and heavy species
(atomic mass > 10).  If any combinations of ions are omitted, then a
warning will be printed once.
The force on the heavy ion is:

.. math::

   \begin{aligned}
   F_z =& \beta \nabla_{||}T_i \\
   \beta =& \frac{3\left(\mu + 5\sqrt{2}Z^2\left(1.1\mu^{5/2} - 0.35\mu^{3/2}\right) - 1\right)}{2.6 - 2\mu + 5.4\mu^2} \\
   \mu =& m_z / \left(m_z + m_i\right)
   \end{aligned}

where subscripts :math:`z` refer to the heavy ion, and :math:`i`
refers to the light ion. The force on the light ion fluid is equal and
opposite: :math:`F_i = -F_z`.

The implementation is in the `ThermalForce` class:

.. doxygenstruct:: ThermalForce
   :members:

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
and then include the group name in the components list

.. code-block:: ini

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

Hydrogen
~~~~~~~~

Multiple isotopes of hydrogen can be evolved, so to keep track of this the
species labels `h`, `d` and `t` are all handled by the same hydrogen atomic
rates calculation. The following might therefore be used

.. code-block:: ini
  
  [hermes]
  components = d, t, reactions

  [reactions]
  type = (
          d + e -> d+ + 2e,  # Deuterium ionisation
          t + e -> t+ + 2e,  # Tritium ionisation
         )

+------------------+---------------------------------------+
| Reaction         | Description                           |
+==================+=======================================+
| h + e -> h+ + 2e | Hydrogen ionisation (Amjuel 2.1.5)    |
+------------------+---------------------------------------+
| d + e -> d+ + 2e | Deuterium ionisation (Amjuel 2.1.5)   |
+------------------+---------------------------------------+
| t + e -> t+ + 2e | Tritium ionisation (Amjuel 2.1.5)     |
+------------------+---------------------------------------+
| h + h+ -> h+ + h | Hydrogen charge exchange              |
+------------------+---------------------------------------+
| d + d+ -> d+ + d | Deuterium charge exchange             |
+------------------+---------------------------------------+
| t + t+ -> t+ + t | Tritium charge exchange               |
+------------------+---------------------------------------+
| h + d+ -> h+ + d | Mixed hydrogen isotope CX             |
+------------------+---------------------------------------+
| d + h+ -> d+ + h |                                       |
+------------------+---------------------------------------+
| h + t+ -> h+ + t |                                       |
+------------------+---------------------------------------+
| t + h+ -> t+ + h |                                       |
+------------------+---------------------------------------+
| d + t+ -> d+ + t |                                       |
+------------------+---------------------------------------+
| t + d+ -> t+ + d |                                       |
+------------------+---------------------------------------+
| h+ + e -> h      | Hydrogen recombination (Amjuel 2.1.8) |
+------------------+---------------------------------------+
| d+ + e -> d      | Deuterium recombination (Amjuel 2.1.8)|
+------------------+---------------------------------------+
| t+ + e -> t      | Tritium recombination (Amjuel 2.1.8)  |
+------------------+---------------------------------------+

The code to calculate the charge exchange rates is in
`hydrogen_charge_exchange.[ch]xx`. This implements reaction 0.1T from
Amjuel (p38), scaled to different isotope masses and finite neutral
particle temperatures by using the effective temperature (Amjuel p43):

.. math::

   T_{eff} = \frac{M}{M_1}T_1 + \frac{M}{M_2}T_2


The effective hydrogenic ionisation rates are calculated using Amjuel
reaction 2.1.5, by D.Reiter, K.Sawada and T.Fujimoto (2016).
Effective recombination rates, which combine radiative and 3-body contributions,
are calculated using Amjuel reaction 2.1.8.

.. doxygenstruct:: HydrogenChargeExchange
   :members:


Helium
~~~~~~

+----------------------+------------------------------------------------------------+
| Reaction             | Description                                                |
+======================+============================================================+
| he + e -> he+ + 2e   | He ionisation, unresolved metastables (Amjuel 2.3.9a)      |
+----------------------+------------------------------------------------------------+
| he+ + e -> he        | He+ recombination, unresolved metastables (Amjuel 2.3.13a) |
+----------------------+------------------------------------------------------------+

The implementation of these rates are in the `AmjuelHeIonisation01`
and `AmjuelHeRecombination10` classes:

.. doxygenstruct:: AmjuelHeIonisation01
   :members:

.. doxygenstruct:: AmjuelHeRecombination10
   :members:

Neon
~~~~

These rates are taken from ADAS (96): SCD and PLT are used for the ionisation
rate and radiation energy loss; ACD and PRB for the recombination rate and radiation
energy loss; and CCD (89) for the charge exchange coupling to hydrogen.
The ionisation potential is also included as a source or sink of energy
for the electrons.

+------------------------+-------------------------------------+
| Reaction               | Description                         |
+========================+=====================================+
| ne + e -> ne+ + 2e     | Neon ionisation                     |
+------------------------+-------------------------------------+
| ne+ + e -> ne+2 + 2e   |                                     |
+------------------------+-------------------------------------+
| ne+2 + e -> ne+3 + 2e  |                                     |
+------------------------+-------------------------------------+
| ne+3 + e -> ne+4 + 2e  |                                     |
+------------------------+-------------------------------------+
| ne+4 + e -> ne+5 + 2e  |                                     |
+------------------------+-------------------------------------+
| ne+5 + e -> ne+6 + 2e  |                                     |
+------------------------+-------------------------------------+
| ne+6 + e -> ne+7 + 2e  |                                     |
+------------------------+-------------------------------------+
| ne+7 + e -> ne+8 + 2e  |                                     |
+------------------------+-------------------------------------+
| ne+8 + e -> ne+9 + 2e  |                                     |
+------------------------+-------------------------------------+
| ne+9 + e -> ne+10 + 2e |                                     |
+------------------------+-------------------------------------+
| ne+ + e -> ne          | Neon recombination                  |
+------------------------+-------------------------------------+
| ne+2 + e -> ne+        |                                     |
+------------------------+-------------------------------------+
| ne+3 + e -> ne+2       |                                     |
+------------------------+-------------------------------------+
| ne+4 + e -> ne+3       |                                     |
+------------------------+-------------------------------------+
| ne+5 + e -> ne+4       |                                     |
+------------------------+-------------------------------------+
| ne+6 + e -> ne+5       |                                     |
+------------------------+-------------------------------------+
| ne+7 + e -> ne+6       |                                     |
+------------------------+-------------------------------------+
| ne+8 + e -> ne+7       |                                     |
+------------------------+-------------------------------------+
| ne+9 + e -> ne+8       |                                     |
+------------------------+-------------------------------------+
| ne+10 + e -> ne+9      |                                     |
+------------------------+-------------------------------------+
| ne+ + h -> ne + h+     | Charge exchange with hydrogen       |
+------------------------+-------------------------------------+
| ne+2 + h -> ne+ + h+   |                                     |
+------------------------+-------------------------------------+
| ne+3 + h -> ne+2 + h+  |                                     |
+------------------------+-------------------------------------+
| ne+4 + h -> ne+3 + h+  |                                     |
+------------------------+-------------------------------------+
| ne+5 + h -> ne+4 + h+  |                                     |
+------------------------+-------------------------------------+
| ne+6 + h -> ne+5 + h+  |                                     |
+------------------------+-------------------------------------+
| ne+7 + h -> ne+6 + h+  |                                     |
+------------------------+-------------------------------------+
| ne+8 + h -> ne+7 + h+  |                                     |
+------------------------+-------------------------------------+
| ne+9 + h -> ne+8 + h+  |                                     |
+------------------------+-------------------------------------+
| ne+10 + h -> ne+9 + h+ |                                     |
+------------------------+-------------------------------------+
| ne+ + d -> ne + d+     | Charge exchange with deuterium      |
+------------------------+-------------------------------------+
| ne+2 + d -> ne+ + d+   |                                     |
+------------------------+-------------------------------------+
| ne+3 + d -> ne+2 + d+  |                                     |
+------------------------+-------------------------------------+
| ne+4 + d -> ne+3 + d+  |                                     |
+------------------------+-------------------------------------+
| ne+5 + d -> ne+4 + d+  |                                     |
+------------------------+-------------------------------------+
| ne+6 + d -> ne+5 + d+  |                                     |
+------------------------+-------------------------------------+
| ne+7 + d -> ne+6 + d+  |                                     |
+------------------------+-------------------------------------+
| ne+8 + d -> ne+7 + d+  |                                     |
+------------------------+-------------------------------------+
| ne+9 + d -> ne+8 + d+  |                                     |
+------------------------+-------------------------------------+
| ne+10 + d -> ne+9 + d+ |                                     |
+------------------------+-------------------------------------+
| ne+ + t -> ne + t+     | Charge exchange with tritium        |
+------------------------+-------------------------------------+
| ne+2 + t -> ne+ + t+   |                                     |
+------------------------+-------------------------------------+
| ne+3 + t -> ne+2 + t+  |                                     |
+------------------------+-------------------------------------+
| ne+4 + t -> ne+3 + t+  |                                     |
+------------------------+-------------------------------------+
| ne+5 + t -> ne+4 + t+  |                                     |
+------------------------+-------------------------------------+
| ne+6 + t -> ne+5 + t+  |                                     |
+------------------------+-------------------------------------+
| ne+7 + t -> ne+6 + t+  |                                     |
+------------------------+-------------------------------------+
| ne+8 + t -> ne+7 + t+  |                                     |
+------------------------+-------------------------------------+
| ne+9 + t -> ne+8 + t+  |                                     |
+------------------------+-------------------------------------+
| ne+10 + t -> ne+9 + t+ |                                     |
+------------------------+-------------------------------------+

The implementation of these rates is in `ADASNeonIonisation`, 
`ADASNeonRecombination` and `ADASNeonCX` template classes:

.. doxygenstruct:: ADASNeonIonisation
   :members:

.. doxygenstruct:: ADASNeonRecombination
   :members:

.. doxygenstruct:: ADASNeonCX
   :members:

Electromagnetic fields
----------------------

These are components which calculate the electric and/or magnetic
fields.

.. _vorticity:

vorticity
~~~~~~~~~

Evolves a vorticity equation, and at each call to transform() uses a matrix
inversion to calculate potential from vorticity.

In this component the Boussinesq approximation is made, so the vorticity equation solved is

.. math::

   \nabla\cdot\left(\frac{\overline{A}\overline{n}}{B^2}\nabla_\perp \phi + \sum_i\frac{A_i}{B^2}\nabla_\perp p_i\right) = \Omega

Where the sum is over species, :math:`\overline{A}` is the average ion
atomic number, and :math:`\overline{n}` is the normalisation density
(i.e. goes to 1 in the normalised equations).  This is a simplified
version of the full expression which is:

.. math::

   \nabla\cdot\left(\sum_i \frac{A_i n_i}{B^2}\nabla_\perp \phi + \sum_i \frac{A_i}{B^2}\nabla_\perp p_i\right) = \Omega

and is derived by replacing

.. math::

   \sum_i A_i n_i \rightarrow \overline{A}\overline{n}

.. doxygenstruct:: Vorticity
   :members:

relax_potential
~~~~~~~~~~~~~~~

This component evolves a vorticity equation, similar to the ``vorticity`` component.
Rather than inverting an elliptic equation at every timestep, this component evolves
the potential in time as a diffusion equation.

.. doxygenstruct:: RelaxPotential
   :members:
