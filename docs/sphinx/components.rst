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
time integration solver. The species charge and atomic mass must be set,
and the initial density should be specified in its own section:

.. code-block:: ini

   [d]
   type = evolve_density, ...

   AA = 2 # Atomic mass
   charge = 0

   [Nd]
   function = 1 - 0.5x # Initial condition, normalised to Nnorm

The equation solved is:

.. math::

   \frac{\partial n}{\partial t} = \nabla\cdot\left[n \left(\frac{1}{B}\mathbf{b}\times\nabla\phi + v_{||}\mathbf{b}\right)\right] + S_n

where the source :math:`S_n` is a combination of external source, and
other processes that nay be included, including drift terms
(e.g. magnetic drift) or atomic processes (e.g. ionization).

Notes:

1. The density will be saved in the output file as `N` + species
   label, e.g `Nd` in the above example.
2. If `diagnose=true` is set in the species options then the net
   source :math:`S_n` is saved as `SN` + species, e.g. `SNd`; the
   external source is saved as `S` + species + `_src` e.g. `Sd_src`.
   The time derivative of density is saved as `ddt(N` + species + `)`
   e.g. `ddt(Nd)`.
3. The density source can be set in the input mesh file as a field
   `S` + species + `_src` e.g. `Sd_src`. This can be overridden by
   specifying the source in the input options.
4. The `poloidal_flows` switch controls whether the X-Y components of
   the ExB flow are included (default is true). If set to `false` then
   ExB flows are only in the X-Z plane.

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
The source, e.g. ``Sd+_feedback``, is calculated as a product of the control signal ``density_source_multiplier``, 
and the array ``density_source_shape`` which defines the source region.
The signal is non-dimensional and the controller depends on the value of ``density_source_shape`` to have a good initial guess of the source.
It should be set to a reasonable value in the units of ``[m-3s-1]``. 
A good reasonable value is the expected steady state domain particle loss (for example due to unrecycled ions at the target).


For example:

.. code-block:: ini

   [d+]
   type = ..., upstream_density_feedback

   density_upstream = 1e19      # Density in m^-3
   density_controller_p = 1e-2  # Feedback controller proportional (p) parameter
   density_controller_i = 1e-3  # Feedback controller integral (i) parameter

   [Nd+]
   source_shape = h(pi - y) * 1e20  # Source shape

There are two additional settings which can make the controller more robust without excessive tuning:

``density_source_positive`` ensures the controller never takes particles away, which can prevent oscillatory
behaviour. Note that this requires some other domain particle sink to ensure control, or else the particle count can never reduce.

``density_integral_positive`` This makes sure the integral component only adds particles. 
The integral component takes a long time to change value, which can result in large overshoots if the initial guess was too small.
This setting mitigates this by disabling the integral term if the density is above the desired value.

Notes:
   - The example cases have their PI parameters tuned properly without the need of the above two settings.
   - Under certain conditions, the use of the PI controller can make the upstream density enter a very small oscillation (~0.05% of upstream value).
   - There is a separate `source` setting that includes a fixed (non varying) density source.

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

Sets the temperature of a species to a fixed value which is constant
in space and time. If the species density is set then this component
also calculates the pressure.

By default only saves the temperature once as a non-evolving variable.
If ``diagnose`` is set then pressure is also saved as a time-evolving
variable.

.. code-block:: ini

   [e]
   type = ..., isothermal

   temperature = 10   # Constant temperature [eV]

.. doxygenstruct:: Isothermal
   :members:


fixed_temperature
~~~~~~~~~~~~~~~~~

Sets the temperature of a species to a fixed value which is constant
in time but can vary in space. If the species density is set then this
component also calculates the pressure.

By default only saves the temperature once as a non-evolving variable.
If ``diagnose`` is set then pressure is also saved as a time-evolving
variable.

.. code-block:: ini

   [e]
   type = ..., fixed_temperature

   temperature = 10 - x   # Spatially dependent temperature [eV]

.. doxygenstruct:: FixedTemperature
   :members:

.. _evolve_pressure:

evolve_pressure
~~~~~~~~~~~~~~~

Evolves the pressure in time. This pressure is named `P<species>` where `<species>`
is the short name of the evolving species e.g. `Pe`.

By default parallel thermal conduction is included, which requires a collision
time. If collisions are not calculated, then thermal conduction should be turned off
by setting `thermal_conduction = false` in the input options.

If the component option ``diagnose = true`` then additional fields
will be saved to the dump files: The species temperature ``T + name``
(e.g. ``Td+`` or ``Te``), the time derivative ``ddt(P + name)``
(e.g. ``ddt(Pd+)`` or ``ddt(Pe)``), and the source of pressure from
other components is saved as ``SP + name`` (e.g. ``SPd+`` or ``SPe``).
The pressure source is the energy density source multiplied by ``2/3``
(i.e. assumes a monatomic species).

.. math::

   \frac{\partial P}{\partial t} = -\nabla\cdot\left(P\mathbf{v}\right) - \frac{2}{3} P \nabla\cdot\mathbf{b}v_{||} + \frac{2}{3}\nabla\cdot\left(\kappa_{||}\mathbf{b}\mathbf{b}\cdot\nabla T\right) + \frac{2}{3}S_E + S_N\frac{1}{2}mNV^2

where :math:`S_E` is the ``energy_source`` (thermal energy source),
and :math:`S_N` is the density source.

Notes:

- Heat conduction through the boundary is turned off currently. This is because
  heat losses are usually calculated at the sheath, so any additional heat conduction
  would be in addition to the sheath heat transmission already included.

The implementation is in `EvolvePressure`:

.. doxygenstruct:: EvolvePressure
   :members:

.. _evolve_energy:

evolve_energy
~~~~~~~~~~~~~

*Note* This is currently under development and has some unresolved
issues with boundary conditions.  Only for testing purposes.

This evolves the sum of species internal energy and parallel kinetic
energy, :math:`\mathcal{E}`:

.. math::

   \mathcal{E} = \frac{1}{\gamma - 1} P + \frac{1}{2}m nv_{||}^2

Note that this component requires the parallel velocity :math:`v_{||}`
to calculate the pressure. It must therefore be listed after a component
that sets the velocity, such as `evolve_momentum`:

.. code-block:: ini

   [d]
   type = ..., evolve_momentum, evolve_energy

The energy density will be saved as `E<species>` (e.g `Ed`) and the
pressure as `P<species>` (e.g. `Pd`). Additional diagnostics, such as the
temperature, can be saved by setting the option `diagnose = true`.

.. doxygenstruct:: EvolveEnergy
   :members:

SNB nonlocal heat flux
~~~~~~~~~~~~~~~~~~~~~~

Calculates the divergence of the electron heat flux using the
Shurtz-Nicolai-Busquet (SNB) model. Uses the BOUT++ implementation which is
`documented here <https://bout-dev.readthedocs.io/en/latest/user_docs/nonlocal.html?#snb-model>`_.

.. doxygenstruct:: SNBConduction
   :members:


Species parallel dynamics
-------------------------

fixed_velocity
~~~~~~~~~~~~~~

Sets the velocity of a species to a fixed value which is constant
in time but can vary in space. If the species density is set then this
component also calculates the momentum.

Saves the temperature once as a non-evolving variable.

.. code-block:: ini

   [e]
   type = ..., fixed_velocity

   velocity = 10 + sin(z)   # Spatially dependent velocity [m/s]

.. doxygenstruct:: FixedVelocity
   :members:


.. _evolve_momentum:

evolve_momentum
~~~~~~~~~~~~~~~

Evolves the momentum `NV<species>` in time. The evolving quantity includes the atomic
mass number, so should be divided by `AA` to obtain the particle flux.

If the component option ``diagnose = true`` then additional fields
will be saved to the dump files: The velocity ``V + name``
(e.g. ``Vd+`` or ``Ve``), the time derivative ``ddt(NV + name)``
(e.g. ``ddt(NVd+)`` or ``ddt(NVe)``), and the source of momentum
density (i.e force density) from other components is saved as ``SNV +
name`` (e.g. ``SNVd+`` or ``SNVe``).

The implementation is in ``EvolveMomentum``:

.. doxygenstruct:: EvolveMomentum
   :members:


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

electron_viscosity
------------------

Calculates the Braginskii electron parallel viscosity, adding a force (momentum source)
to the electron momentum equation:

.. math::

   F = \sqrt{B}\nabla\cdot\left[\frac{\eta_e}{B}\mathbf{b}\mathbf{b}\cdot\nabla\left(\sqrt{B}V_{||e}\right)\right]

The electron parallel viscosity is

.. math::

   \eta_e = \frac{4}{3} 0.73 p_e \tau_e

where :math:`\tau_e` is the electron collision time. The collisions between electrons
and all other species therefore need to be calculated before this component is run:

.. code-block:: ini

   [hermes]
   components = ..., e, ..., collisions, electron_viscosity

.. doxygenstruct:: ElectronViscosity
   :members:

ion_viscosity
-------------

Adds ion viscosity terms to all charged species that are not electrons.
The collision frequency is required so this is a top-level component that
must be calculated after collisions:

.. code-block:: ini

   [hermes]
   components =  ..., collisions, ion_viscosity

By default only the parallel diffusion of momentum is included, adding a force to each
ion's momentum equation:

.. math::

   F = \sqrt{B}\nabla\cdot\left[\frac{\eta_i}{B}\mathbf{b}\mathbf{b}\cdot\nabla\left(\sqrt{B}V_{||i}\right)\right]

The ion parallel viscosity is

.. math::

   \eta_i = \frac{4}{3} 0.96 p_i \tau_i

If the `perpendicular` option is set:

.. code-block:: ini

   [ion_viscosity]
   perpendicular = true # Include perpendicular flows

Then the ion scalar viscous pressure is calculated as:

.. math::

   \Pi_{ci} = \Pi_{ci||} + \Pi_{ci\perp}

where :math:`\Pi_{ci||}` corresponds to the parallel diffusion of momentum above.

.. math::

   \Pi_{ci||} = - 0.96 \frac{2p_i\tau_i}{\sqrt{B}} \partial_{||}\left(\sqrt{B} V_{||i}\right)

The perpendicular part is calculated from:

.. math::

   \begin{aligned}\Pi_{ci\perp} =& 0.96 p_i\tau_i \kappa \cdot \left[\mathbf{V}_E + \mathbf{V}_{di} + 1.16\frac{\mathbf{b}\times\nabla T_i}{B} \right] \\
   =& -0.96 p_i\tau_i\frac{1}{B}\left(\mathbf{b}\times\kappa\right)\cdot\left[\nabla\phi + \frac{\nabla p_i}{en_i} + 1.61\nabla T_i \right]\end{aligned}


A parallel force term is added, in addition to the parallel viscosity above:

.. math::

   F = -\frac{2}{3}B^{3/2}\partial_{||}\left(\frac{\Pi_{ci\perp}}{B^{3/2}}\right)
   
In the vorticity equation the viscosity appears as a divergence of a current:

.. math::

   \mathbf{J}_{ci} = \frac{\Pi_{ci}}{2}\nabla\times\frac{\mathbf{b}}{B} - \frac{1}{3}\frac{\mathbf{b}\times\nabla\Pi_{ci}}{B}

that transfers energy between ion internal energy and :math:`E\times B` energy:

.. math::

   \begin{aligned}\frac{\partial \omega}{\partial t} =& \ldots + \nabla\cdot\mathbf{J}_{ci} \\
   \frac{\partial p_i}{\partial t} =& \ldots - \mathbf{J}_{ci}\cdot\nabla\left(\phi + \frac{p_i}{n_0}\right)\end{aligned}

Note that the sum of the perpendicular and parallel contributions to the ion viscosity act to damp
the net poloidal flow. This can be seen by assuming that :math:`\phi`, :math:`p_i` and :math:`T_i`
are flux functions. We can then write:

.. math::

   \Pi_{ci\perp} = -0.96 p_i\tau_i \frac{1}{B}\left(\mathbf{b}\times\kappa\right)\cdot\nabla\psi F\left(\psi\right)

where

.. math::

   F\left(\psi\right) = \frac{\partial\phi}{\partial\psi} + \frac{1}{en}\frac{\partial p_i}{\partial\psi} + 1.61\frac{\partial T_i}{\partial\psi}

Using the approximation

.. math::

   \left(\mathbf{b}\times\kappa\right)\cdot\nabla\psi \simeq -RB_\zeta \partial_{||}\ln B

expanding:

.. math::

   \frac{2}{\sqrt{B}}\partial_{||}\left(\sqrt{B}V_{||i}\right) = 2\partial_{||}V_{||i} + V_{||i}\partial_{||}\ln B

and neglecting parallel gradients of velocity gives:

.. math::

   \Pi_{ci} \simeq 0.96 p_i\tau_i \left[ \frac{RB_{\zeta}}{B}F\left(\psi\right) - V_{||i} \right]\partial_{||}\ln B

   
**Notes** and implementation details:
- The magnitude of :math:`\Pi_{ci\perp}` and :math:`\Pi_{ci||}` are individually
  limited to be less than or equal to the scalar pressure :math:`Pi` (though can have
  opposite sign). The reasoning is that if these off-diagonal terms become large then
  the model is likely breaking down. Occasionally happens in low-density regions.

   
.. doxygenstruct:: IonViscosity
   :members:

simple_conduction
-----------------

This is a simplified parallel heat conduction model that can be used when a linearised model is needed.
If used, the thermal conduction term in `evolve_pressure` component should be disabled.

.. code-block:: ini

   [hermes]
   components = e, ...

   [e]
   type = evolve_pressure, simple_conduction

   thermal_conduction = false  # Disable term in evolve_pressure

To linearise the heat conduction the temperature and density used in
calculating the Coulomb logarithm and heat conduction coefficient can
be fixed by specifying `conduction_temperature` and
`conduction_density`.

Note: For hydrogenic plasmas this produces very similar parallel electron
heat conduction as the `evolve_pressure` term with electron-electron collisions
disabled.

.. doxygenstruct:: SimpleConduction
   :members:

Drifts
------

The ExB drift is included in the density, momentum and pressure evolution equations if
potential is calculated. Other drifts can be added with the following components.

diamagnetic_drift
~~~~~~~~~~~~~~~~~

Adds diamagnetic drift terms to all species' density, pressure and parallel momentum
equations. Calculates the diamagnetic drift velocity as

.. math::

   \mathbf{v}_{dia} = \frac{T}{q} \nabla\times\left(\frac{\mathbf{b}}{B}\right)

where the curvature vector :math:`\nabla\times\left(\frac{\mathbf{b}}{B}\right)`
is read from the `bxcv` mesh input variable.

.. doxygenstruct:: DiamagneticDrift
   :members:


polarisation_drift
~~~~~~~~~~~~~~~~~~

This calculates the polarisation drift of all charged species,
including ions and electrons. It works by approximating the drift
as a potential flow:

.. math::

   \mathbf{v}_{pol} = - \frac{m}{q B^2} \nabla_\perp\phi_{pol}

where :math:`\phi_{pol}` is approximately the time derivative of the
electrostatic potential :math:`\phi` in the frame of the fluid, with
an ion diamagnetic contribution. This is calculated by inverting a
Laplacian equation similar to that solved in the vorticity equation.

This component needs to be run after all other currents have been
calculated.  It marks currents as used, so out-of-order modifications
should raise errors.

See the `examples/blob2d-vpol` example, which contains:

.. code-block:: ini

   [hermes]
   components = e, vorticity, sheath_closure, polarisation_drift

   [polarisation_drift]
   diagnose = true

Setting `diagnose = true` saves `DivJ` to the dump files with the divergence of all
currents except polarisation, and `phi_pol` which is the polarisation flow potential.

.. doxygenstruct:: PolarisationDrift
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

Boundary conditions
-------------------

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

.. _neutral_boundary:

neutral_boundary
~~~~~~~~~~~~~~~~

Sets Y (sheath/target) boundary conditions on neutral particle
density, temperature and pressure. A no-flow boundary condition
is set on parallel velocity and momentum. It is a species-specific
component and so goes in the list of components for the species
that the boundary condition should be applied to.

An energy sink is added to the flux of heat to the wall, with
heat flux `q`:

.. math::

   q = \gamma_{heat} n T v_{th}

   v_{th} = \sqrt{eT / m}

The factor `gamma_heat`

.. code-block:: ini

   [hermes]
   components = d

   [d]
   type = ... , neutral_boundary

   gamma_heat = 3  # Neutral boundary heat transmission coefficient
   neutral_lower_y = true  # Boundary on lower y?
   neutral_upper_y = true  # Boundary on upper y?

.. doxygenstruct:: NeutralBoundary
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

.. code-block:: ini

   [hermes]
   components = ... , collisions, neutral_parallel_diffusion

   [neutral_parallel_diffusion]
   dneut = 1         # Diffusion multiplication factor
   diagnose = true   # This enables diagnostic output for each species


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

For collisions between charged particles. In the following all
quantities are in SI units except the temperatures: :math:`T` is in
eV, so :math:`eT` has units of Joules.

Debye length :math:`\lambda_D`

.. math::

   \lambda_D = \sqrt{\frac{\epsilon_0 T_e}{n_e e}}
   
Coulomb logarithm, from [NRL formulary 2019], adapted to SI units

- For thermal electron-electron collisions

  .. math::

     \ln \lambda_{ee} = 30.4 - \frac{1}{2} \ln\left(n_e\right) + \frac{5}{4}\ln\left(T_e\right) - \sqrt{10^{-5} + \left(\ln T_e - 2\right)^2 / 16} 

  where the coefficient (30.4) differs from the NRL value due to
  converting density from cgs to SI units (:math:`30.4 = 23.5 -
  0.5\ln\left(10^{-6}\right)`).


- Electron-ion collisions

  .. math::

     \ln \lambda_{ei} = \left\{\begin{array}{ll}
                              10 & \textrm{if } T_e < 0.1 \textrm{eV or } n_e < 10^{10}m^{-3} \\
                              30 - \frac{1}{2}\ln\left(n_e\right) - \ln(Z) + \frac{3}{2}\ln\left(T_e\right) & \textrm{if } T_im_e/m_i < T_e < 10Z^2 \\
                              31 - \frac{1}{2}\ln\left(n_e\right) + \ln\left(T_e\right) & \textrm{if } T_im_e/m_i < 10Z^2 < T_e \\
                              23 - \frac{1}{2}\ln\left(n_i\right) + \frac{3}{2}\ln\left(T_i\right) - \ln\left(Z^2\mu\right) & \textrm{if } T_e < T_im_e/m_i \\
                              \end{array}\right.
     
- Mixed ion-ion collisions
  
  .. math::

     \ln \lambda_{ii'} = 29.91 - ln\left[\frac{ZZ'\left(\mu + \mu'\right)}{\mu T_{i'} + \mu'T_i}\left(\frac{n_iZ^2}{T_i} + \frac{n_{i'} Z'^2}{T_{i'}}\right)^{1/2}\right]

  where like the other expressions the different constant is due to
  converting from cgs to SI units: :math:`29.91 = 23 -
  0.5\ln\left(10^{-6}\right)`.

The frequency of charged species `a` colliding with charged species `b` is

.. math::

   \nu_{ab} = \frac{1}{3\pi^{3/2}\epsilon_0^2}\frac{Z_a^2 Z_b^2 n_b \ln\Lambda}{\left(v_a^2 + v_b^2\right)^{3/2}}\frac{\left(1 + m_a / m_b\right)}{m_a^2}


Note that the cgs expression in Hinton is divided by :math:`\left(4\pi\epsilon_0\right)^2` to get
the expression in SI units. The thermal speeds in this expression are defined as:

.. math::

   v_a^2 = 2 e T_a / m_a

Note that with this definition we recover the `Braginskii expressions
<https://farside.ph.utexas.edu/teaching/plasma/lectures1/node35.html>`_
for e-i and i-i collision times.

For conservation of momentum, the collision frequencies :math:`\nu_{ab}` and :math:`\nu_{ba}` are
related by:

.. math::

   m_a n_a \nu_{ab} = m_b n_b \nu_{ba}

Momentum exchange, force on species `a` due to collisions with species `b`:

.. math::

   F_{ab} = C_m \nu_{ab} m_a n_a \left( u_b - u_a \right)

Where the coefficient :math:`C_m` for parallel flows depends on the species: For most combinations
of species this is set to 1, but for electron-ion collisions the Braginskii coefficients are used:
:math:`C_m = 0.51` if ion charge :math:`Z_i = 1`;  0.44 for :math:`Z_i = 2`; 0.40 for :math:`Z_i = 3`;
and 0.38 is used for :math:`Z_i \ge 4`. Note that this coefficient should decline further with
increasing ion charge, tending to 0.29 as :math:`Z_i \rightarrow \infty`.

Frictional heating is included by default, but can be disabled by
setting the `frictional_heating` option to `false`. When enabled it
adds a source of thermal energy corresponding to the resistive heating
term:

.. math::

   Q_{ab,F} = \frac{m_b}{m_a + m_b} \left( u_b - u_a \right) F_{ab}

This term has some important properties:

1. It is always positive: Collisions of two species with the same
   temperature never leads to cooling.
2. It is Galilean invariant: Shifting both species' velocity by the
   same amount leaves :math:`Q_{ab,F}` unchanged.
3. If both species have the same mass, the thermal energy
   change due to slowing down is shared equally between them.
4. If one species is much heavier than the other, for example
   electron-ion collisions, the lighter species is preferentially
   heated. This recovers e.g. Braginskii expressions for :math:`Q_{ei}`
   and :math:`Q_{ie}`.

This can be derived by considering the exchange of energy
:math:`W_{ab,F}` between two species at the same temperature but
different velocities. If the pressure is evolved then it contains
a term that balances the change in kinetic energy due to changes
in velocity:

.. math::

   \begin{aligned}
   \frac{\partial}{\partial t}\left(m_a n_a u_a\right) =& \ldots + F_{ab} \\
   \frac{\partial}{\partial t}\left(\frac{3}{2}p_a\right) =& \ldots - F_{ab} u_a + W_{ab, F}
   \end{aligned}

For momentum and energy conservation we must have :math:`F_{ab}=-F_{ba}`
and :math:`W_{ab,F} = -W_{ba,F}`. Comparing the above to the
`Braginskii expression
<https://farside.ph.utexas.edu/teaching/plasma/lectures/node35.html>`_
we see that for ion-electron collisions the term :math:`- F_{ab}u_a + W_{ab, F}`
goes to zero, so :math:`W_{ab, F} \sim u_aF_{ab}` for
:math:`m_a \gg m_b`. An expression that has all these desired properties
is

.. math::

   W_{ab,F} = \left(\frac{m_a u_a + m_b u_a}{m_a + m_b}\right)F_{ab}

which is not Galilean invariant but when combined with the :math:`- F_{ab} u_a`
term gives a change in pressure that is invariant, as required.
   
Thermal energy exchange, heat transferred to species :math:`a` from
species :math:`b` due to temperature differences, is given by:

.. math::

   Q_{ab,T} = \nu_{ab}\frac{3n_a m_a\left(T_b - T_a\right)}{m_a + m_b}

- Ion-neutral and electron-neutral collisions

  The cross-section for elastic collisions between charged and neutral
  particles can vary significantly. Here for simplicity we just take
  a value of :math:`5\times 10^{-19}m^2` from the NRL formulary.

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

.. _recycling:

recycling
~~~~~~~~~

This component calculates the flux of a species into a Y boundary,
due to recycling of flow out of the boundary of another species.

The boundary fluxes might be set by sheath boundary conditions,
which potentially depend on the density and temperature of all species.
Recycling therefore can't be calculated until all species boundary conditions
have been set. It is therefore expected that this component is a top-level
component which comes after boundary conditions are set.

The recycling component has a `species` option, that is a list of species
to recycle. For each of the species in that list, `recycling` will look in
the corresponding section for the options `recycle_as`, `recycle_multiplier`
and `recycle_energy`.

For example, recycling `d+` ions into `d` atoms with a recycling fraction
of 1. Each returning atom has an energy of 3.5eV:

.. code-block:: ini

   [hermes]
   components = d+, d, sheath_boundary, recycling

   [recycling]
   species = d+   # Comma-separated list of species to recycle

   [d+]
   recycle_as = d         # Species to recycle as
   recycle_multiplier = 1 # Recycling fraction
   recycle_energy = 3.5   # Energy of recycled particles [eV]

.. doxygenstruct:: Recycling
   :members:

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

Notes:

1. Charge exchange channel diagnostics: For two species `a` and `b`,
   the channel `Fab_cx` is a source of momentum for species `a` due to
   charge exchange with species `b`. There are corresponding sinks for
   the products of the charge exchange reaction which are not saved.

   For example,reaction `d + t+ -> d+ + t` will save the following
   forces (momentum sources):
   - `Fdt+_cx` is a source of momentum for deuterium atoms `d` and sink of momentum for deuterium ions `d+`.
   - `Ft+d_cx` is a source of momentum for tritium ions `t+` and sink of momentum for tritium atoms `t`

   The reason for this convention is the existence of the inverse reactions:
   `t + d+ -> t+ + d` outputs diagnostics `Ftd+_cx` and `Fd+t_cx`.

2. Reactions typically convert species from one to another, leading to
   a transfer of mass momentum and energy. For a reaction converting
   species :math:`a` to species :math:`b` at rate :math:`R` (units
   of events per second per volume) we have transfers:

   .. math::

      \begin{aligned}
      \frac{\partial}{\partial t} n_a =& \ldots - R \\
      \frac{\partial}{\partial t} n_b =& \ldots + R \\
      \frac{\partial}{\partial t}\left( m n_a u_a\right) =& \ldots + F_{ab} \\
      \frac{\partial}{\partial t}\left( m n_a u_a\right) =& \ldots + F_{ba} \\
      \frac{\partial}{\partial t}\left( \frac{3}{2} p_a \right) =& \ldots - F_{ab}u_a + W_{ab} - \frac{1}{2}mRu_a^2 \\
      \frac{\partial}{\partial t}\left( \frac{3}{2} p_b \right) =& \ldots - F_{ba}u_b + W_{ba} + \frac{1}{2}mRu_b^2
      \end{aligned}
      
  where both species have the same mass: :math:`m_a = m_b = m`. In the
  pressure equations the :math:`-F_{ab}u_a` comes from splitting the
  kinetic and thermal energies; :math:`W_{ab}=-W_{ba}` is the energy
  transfer term that we need to find; The final term balances the loss
  of kinetic energy at fixed momentum due to a particle source or
  sink.

  The momentum transfer :math:`F_{ab}=-F{ba}` is the momentum carried
  by the converted ions: :math:`F_{ab}=-m R u_a`. To find
  :math:`W_{ab}` we note that for :math:`p_a = 0` the change in pressure
  must go to zero: :math:`-F_{ab}u_a + W_{ab} -\frac{1}{2}mRu_a^2 = 0`.

  .. math::

      \begin{aligned}
      W_{ab} =& F_{ab}u_a + \frac{1}{2}mRu_a^2 \\
      =& - mR u_a^2 + \frac{1}{2}mRu_a^2\\
      =& -\frac{1}{2}mRu_a^2
      \end{aligned}

  Substituting into the above gives:

  .. math::

     \begin{aligned}
     \frac{\partial}{\partial t}\left( \frac{3}{2} p_b \right) =& \ldots - F_{ba}u_b + W_{ba} + \frac{1}{2}mRu_b^2 \\
     =& \ldots - mRu_au_b + \frac{1}{2}mRu_a^2 + \frac{1}{2}mRu_a^2 \\
     =& \ldots + \frac{1}{2}mR\left(u_a - u_b\right)^2
     \end{aligned}

  This has the property that the change in pressure of both species is
  Galilean invariant. This transfer term is included in the Amjuel reactions
  and hydrogen charge exchange.
     
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
`hydrogen_charge_exchange.[ch]xx`. This implements reaction 3.1.8 from
Amjuel (p43), scaled to different isotope masses and finite neutral
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

Fixed fraction radiation
~~~~~~~~~~~~~~~~~~~~~~~~

These components produce volumetric electron energy losses, but don't
otherwise modify the plasma solution: Their charge and mass density
are not calculated, and there are no interactions with other species
or boundary conditions.

The ``fixed_fraction_hutchinson_carbon`` component calculates radiation due to carbon
in coronal equilibrium, using a simple formula from `I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994) <https://doi.org/10.1088/0029-5515/34/10/I04>`_:

.. math::

   L\left(T_e\right) = 2\times 10^{-31} \frac{\left(T_e/10\right)^3}{1 + \left(T_e / 10\right)^{4.5}}

which has units of :math:`Wm^3` with :math:`T_e` in eV.

To use this component you can just add it to the list of components and then
configure the impurity fraction:

.. code-block:: ini

   [hermes]
   components = ..., fixed_fraction_hutchinson_carbon, ...

   [fixed_fraction_hutchinson_carbon]
   fraction = 0.05   # 5% of electron density
   diagnose = true   # Saves Rfixed_fraction_carbon to output

Or to customise the name of the radiation output diagnostic a section can be
defined like this:

.. code-block:: ini

   [hermes]
   components = ..., c, ...

   [c]
   type = fixed_fraction_hutchinson_carbon
   fraction = 0.05   # 5% of electron density
   diagnose = true   # Saves Rc (R + section name)


Carbon is also provided as an ADAS rate along with nitrogen, neon and argon. The component names are  
``fixed_fraction_carbon``, ``fixed_fraction_nitrogen``, ``fixed_fraction_neon`` and ``fixed_fraction_argon``.

These can be used in the same way as ``fixed_fraction_hutchinson_carbon``. Each rate is in the form of a 10 coefficient 
log-log polynomial fit of data obtained using the open source tool `radas <https://github.com/cfs-energy/radas>`_.
The :math:`n {\tau}` parameter representing the density and residence time assumed in the radas 
collisional-radiative model has been set to :math:`1\times 10^{20} \times 0.5ms` based on `David Moulton et al 2017 Plasma Phys. Control. Fusion 59(6) <https://doi.org10.1088/1361-6587/aa6b13>`_.

Each rate has an upper and lower bound beyond which the rate remains constant. 
Please refer to the source code in `fixed_fraction_radiation.hxx` for the coefficients and bounds used for each rate.

Electromagnetic fields
----------------------

These are components which calculate the electric and/or magnetic
fields.

.. _vorticity:

vorticity
~~~~~~~~~

Evolves a vorticity equation, and at each call to transform() uses a matrix
inversion to calculate potential from vorticity.

In this component the Boussinesq approximation is made, so the
vorticity equation solved is

.. math::

   \nabla\cdot\left(\frac{\overline{A}\overline{n}}{B^2}\nabla_\perp \phi\right) \underbrace{+ \nabla\cdot\left(\sum_i\frac{A_i}{Z_i B^2}\nabla_\perp p_i\right)}_{\mathrm{if diamagnetic\_polarisation}} = \Omega

Where the sum is over species, :math:`\overline{A}` is the average ion
atomic number, and :math:`\overline{n}` is the normalisation density
(i.e. goes to 1 in the normalised equations). The ion diamagnetic flow
terms in this Boussinesq approximation can be written in terms of an
effective ion pressure :math:`\hat{p}`:

.. math::

   \hat{p} \equiv \sum_i \frac{A_i}{\overline{A} Z_i} p_i

as

.. math::

   \nabla\cdot\left[\frac{\overline{A}\overline{n}}{B^2}\nabla_\perp \left(\phi + \frac{\hat{p}}{\overline{n}}\right) \right] = \Omega
   
Note that if ``diamagnetic_polarisation = false`` then the ion
pressure terms are removed from the vorticity, and also from other ion
pressure terms coming from the polarisation current
(i.e. :math:`\hat{p}\rightarrow 0`.

This is a simplified version of the full vorticity definition which is:

.. math::

   \nabla\cdot\left(\sum_i \frac{A_i n_i}{B^2}\nabla_\perp \phi + \sum_i \frac{A_i}{Z_i B^2}\nabla_\perp p_i\right) = \Omega

and is derived by replacing

.. math::

   \sum_i A_i n_i \rightarrow \overline{A}\overline{n}

In the case of multiple species, this Boussinesq approximation means that the ion diamagnetic flow
terms 

The vorticity equation that is integrated in time is

.. math::

   \begin{aligned}\frac{\partial \Omega}{\partial t} =& \nabla\cdot\left(\mathbf{b}\sum_s Z_s n_sV_{||s}\right) \\
   &+ \underbrace{\nabla\cdot\left(\nabla\times\frac{\mathbf{b}}{B}\sum_s p_s\right)}_{\textrm{if diamagnetic}} + \underbrace{\nabla\cdot\mathbf{J_{exb}}}_{\mathrm{if exb\_advection}} \\
   &+ \nabla\cdot\left(\mathbf{b}J_{extra}\right)\end{aligned}

The nonlinearity :math:`\nabla\cdot\mathbf{J_{exb}}` is part of the
divergence of polarisation current. In its simplified form when
``exb_advection_simplified = true``, this is the :math:`E\times B`
advection of vorticity:

.. math::

   \nabla\cdot\mathbf{J_{exb}} = -\nabla\cdot\left(\Omega \mathbf{V}_{E\times B}\right)

When ``exb_advection_simplified = false`` then the more complete
(Boussinesq approximation) form is used:

.. math::

   \nabla\cdot\mathbf{J_{exb}} = -\nabla\cdot\left[\frac{\overline{A}}{2B^2}\nabla_\perp\left(\mathbf{V}_{E\times B}\cdot\nabla \hat{p}\right) + \frac{\Omega}{2} \mathbf{V}_{E\times B} + \frac{\overline{A}\overline{n}}{2B^2}\nabla_\perp^2\phi\left(\mathbf{V}_{E\times B} + \frac{\mathbf{b}}{B}\times\nabla\hat{p}\right) \right]
   
The form of the vorticity equation is based on `Simakov & Catto
<https://doi.org/10.1063/1.1623492>`_ (corrected in `erratum 2004
<https://doi.org/10.1063/1.1703527>`_), in the Boussinesq limit and
with the first term modified to conserve energy. In the limit of zero
ion pressure and constant :math:`B` it reduces to the simplified form.

.. doxygenstruct:: Vorticity
   :members:

relax_potential
~~~~~~~~~~~~~~~

This component evolves a vorticity equation, similar to the ``vorticity`` component.
Rather than inverting an elliptic equation at every timestep, this component evolves
the potential in time as a diffusion equation.

.. doxygenstruct:: RelaxPotential
   :members:

electromagnetic
~~~~~~~~~~~~~~~

This component modifies the definition of momentum of all species, to
include the contribution from the electromagnetic potential
:math:`A_{||}`.

Assumes that "momentum" :math:`p_s` calculated for all species
:math:`s` is

.. math::

   p_s = m_s n_s v_{||s} + Z_s e n_s A_{||}

which arises once the electromagnetic contribution to the force on
each species is included in the momentum equation. This is normalised
so that in dimensionless quantities

.. math::

   p_s = A n v_{||} + Z n A_{||}

where :math:`A` and :math:`Z` are the atomic number and charge of the
species.

The current density :math:`j_{||}` in SI units is

.. math::

   j_{||} = -\frac{1}{\mu_0}\nabla_\perp^2 A_{||}

which when normalised in Bohm units becomes

.. math::

   j_{||} = - \frac{1}{\beta_{em}}\nabla_\perp^2 A_{||}

where :math:`\beta_{em}` is a normalisation parameter which is half
the plasma electron beta as normally defined:

.. math::

   \beta_{em} = \frac{\mu_0 e \overline{n} \overline{T}}{\overline{B}^2}

To convert the species momenta into a current, we take the sum of
:math:`p_s Z_s e / m_s`. In terms of normalised quantities this gives:

.. math::

   - \frac{1}{\beta_{em}} \nabla_\perp^2 A_{||} + \sum_s \frac{Z^2 n_s}{A}A_{||} = \sum_s \frac{Z}{A} p_s

.. doxygenstruct:: Electromagnetic
   :members:
