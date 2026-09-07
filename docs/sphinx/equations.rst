.. _sec-equations:

Equations
==========

This section contains components which determine which equations are
being solved in the code. There are two broad classes of components:

Whole equations
   For example, ``fixed_temperature``, ``evolve_pressure``,
   ``evolve_energy`` allow the solution of energy in three levels
   of fidelity: constant temperature, a pressure equation and the
   conservative total energy equation. ``neutral_mixed`` contains
   both parallel and perpendicular transport of neutrals and
   has several equations included within.

Terms
   For example, ``anomalous_diffusion`` adds cross-field transport
   to the density, energy and momentum equations if they are available
   while  ``diamagnetic_drift`` and ``polarisation_drift``
   add drift terms.

Please refer to the examples for common component configurations.


Species density
---------------

The density of a species can be calculated in several different ways,
and are usually needed by other components.

.. _sec-fixed_density:

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

.. _sec-evolve_density:

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

   \frac{\partial n}{\partial t} = -\nabla\cdot\left[n \left(\frac{1}{B}\mathbf{b}\times\nabla\phi + v_{||}\mathbf{b}\right)\right] + S_n

where the source :math:`S_n` is a combination of external source, and
other processes that nay be included, including drift terms
(e.g. magnetic drift) or atomic processes (e.g. ionization).

Notes:

1. The density will be saved in the output file as ``N`` + species
   label, e.g ``Nd`` in the above example.
2. If ``diagnose=true`` is set in the species options then the net
   source :math:`S_n` is saved as ``SN`` + species, e.g. ``SNd``; the
   external source is saved as ``S`` + species + ``_src`` e.g. ``Sd_src``.
   The time derivative of density is saved as ``ddt(N`` + species + ``)``
   e.g. ``ddt(Nd)``.
3. The density source can be set in the input mesh file as a field
   ``S`` + species + ``_src`` e.g. ``Sd_src``. This can be overridden by
   specifying the source in the input options.
4. The ``poloidal_flows`` switch controls whether the X-Y components of
   the ExB flow are included (default is true). If set to ``false`` then
   ExB flows are only in the X-Z plane.

The implementation is in the `EvolveDensity` class:

.. doxygenstruct:: EvolveDensity
   :members:



fixed_fraction_ions
~~~~~~~~~~~~~~~~~~~

This sets the density of a species to a fraction of the electron density.

.. doxygenstruct:: FixedFractionIons
   :members:

.. _sec-quasineutral:

quasineutral
~~~~~~~~~~~~

This component sets the density of one species, so that the overall
charge density is zero everywhere. This must therefore be done after
all other charged species densities have been calculated. It only
makes sense to use this component for species with a non-zero charge.

.. doxygenstruct:: Quasineutral
   :members:

Species pressure and temperature
--------------------------------

.. _sec-isothermal:

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

.. _sec-evolve_pressure:

evolve_pressure
~~~~~~~~~~~~~~~

Evolves the pressure in time. This pressure is named ``P<species>`` where ``<species>``
is the short name of the evolving species e.g. ``Pe``.

Parallel conduction is included if the global
:ref:`sec-braginskii_conduction` component has been used.

If the component option ``diagnose = true`` then additional fields
will be saved to the dump files: The species temperature ``T + name``
(e.g. ``Td+`` or ``Te``), the time derivative ``ddt(P + name)``
(e.g. ``ddt(Pd+)`` or ``ddt(Pe)``), and the source of pressure from
other components is saved as ``SP + name`` (e.g. ``SPd+`` or ``SPe``).
The pressure source is the energy density source multiplied by ``2/3``
(i.e. assumes a monatomic species).

.. math::

   \frac{\partial P}{\partial t} = -\nabla\cdot\left(P\mathbf{v}\right) - \frac{2}{3} P \nabla\cdot\mathbf{b}v_{||} + \frac{2}{3}S_E + \frac{2}{3}S_N\frac{1}{2}mNV^2

where :math:`S_E` is the ``energy_source`` (thermal energy source),
and :math:`S_N` is the density source. If conduction has been
calculated, it will be included in the ``energy_source`` term.

Notes:

- Heat conduction through the boundary is turned off currently. This is because
  heat losses are usually calculated at the sheath, so any additional heat conduction
  would be in addition to the sheath heat transmission already included.

The implementation is in `EvolvePressure`:

.. doxygenstruct:: EvolvePressure
   :members:

.. _sec-evolve_energy:

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

The energy density will be saved as ``E<species>`` (e.g ``Ed``) and the
pressure as ``P<species>`` (e.g. ``Pd``). Additional diagnostics, such as the
temperature, can be saved by setting the option ``diagnose = true``.

.. doxygenstruct:: EvolveEnergy
   :members:

SNB nonlocal heat flux
~~~~~~~~~~~~~~~~~~~~~~

Calculates the divergence of the electron heat flux using the
Shurtz-Nicolai-Busquet (SNB) model. Uses the BOUT++ implementation which is
`documented here <https://bout-dev.readthedocs.io/en/latest/user_docs/nonlocal.html?#snb-model>`__.

.. doxygenstruct:: SNBConduction
   :members:


simple_conduction
~~~~~~~~~~~~~~~~~~~~~~

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

.. _sec-braginskii_conduction:

braginskii_conduction
~~~~~~~~~~~~~~~~~~~~~

This is a global component that calculates the parallel thermal
conduction for all species that use :ref:`sec-evolve_pressure` or
:ref:`sec-evolve_energy`, storing it in `energy_source`. If this is not
desired for a particular species then it can be turned off by setting
:code:`thermal_conduction = false` in the input options for that species.

This component requires a collision time to have been calculated
(i.e., with the :ref:`sec-collisions` component). It is
recommended that this be one of the last component to run, to ensure density,
pressure, and temperature have their final values. However, it must be
run before :ref:`sec-recycling`, as that component will need to use the
`energy_flow_ylow` value, to which conduction contributes.

The choice of collision frequency used for conduction is set by the
flag `conduction_collisions_mode`: `multispecies` uses all available
collision frequencies involving the chosen species, while `braginskii`
uses only self-collisions .The default is `multispecies` and it is
recommended for use if solving more than one ion.  If you are solving
for a single ion and want to recover Braginskii, use the `braginskii`
mode.

.. math::
    \nabla\cdot\left(\kappa_{||}\mathbf{b}\mathbf{b}\cdot\nabla T\right)

.. doxygenstruct:: BraginskiiConduction
   :members:


Species parallel dynamics
-------------------------

fixed_velocity
~~~~~~~~~~~~~~

Sets the velocity of a species to a fixed value which is constant
in time but can vary in space. If the species density is set then this
component also calculates the momentum.

Saves the temperature once as a non-evolving variable.

The velocity may be set in the mesh file as an array (2D or 3D), or in
the options as an expression. The options value overrides the value in
the mesh. If neither mesh array nor option are set then an exception
will be thrown. Both mesh array and option should be specified in
units of meters per second.

The name of the array in the mesh file is ``V<name>0`` where
``<name>`` is the name of the species e.g. for species ``e``
(electrons), ``fixed_velocity`` will try to read ``Ve0`` from the mesh
file, and then the ``velocity`` option in the ``[e]`` section:

.. code-block:: ini

   [e]
   type = ..., fixed_velocity

   velocity = 10 + sin(z)   # Spatially dependent velocity [m/s]

.. doxygenstruct:: FixedVelocity
   :members:


.. _sec-evolve_momentum:

evolve_momentum
~~~~~~~~~~~~~~~

Evolves the momentum `NV<species>` in time. The evolving quantity includes the atomic
mass number, so should be divided by `AA` to obtain the particle flux.

If the component option :code:`diagnose = true` then additional fields
will be saved to the dump files: The velocity ``V + name``
(e.g. ``Vd+`` or ``Ve``), the time derivative ``ddt(NV + name)``
(e.g. ``ddt(NVd+)`` or ``ddt(NVe)``), and the source of momentum
density (i.e force density) from other components is saved as ``SNV +
name`` (e.g. ``SNVd+`` or ``SNVe``).

The implementation is in ``EvolveMomentum``:

.. doxygenstruct:: EvolveMomentum
   :members:


.. _sec-zero_current:

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

.. doxygenstruct:: ZeroCurrent
   :members:

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

.. _sec-braginskii_electron_viscosity:

braginskii_electron_viscosity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
   components = ..., e, ..., braginskii_collisions, braginskii_electron_viscosity

.. doxygenstruct:: BraginskiiElectronViscosity
   :members:

.. _sec-braginskii_ion_viscosity:

braginskii_ion_viscosity
~~~~~~~~~~~~~~~~~~~~~~~~

Adds ion viscosity terms to all charged species that are not electrons.
The collision frequency is required so this is a top-level component that
must be calculated after collisions:

.. code-block:: ini

   [hermes]
   components =  ..., braginskii_collisions, braginskii_ion_viscosity

By default only the parallel diffusion of momentum is included, adding a force to each
ion's momentum equation:

.. math::

   F = \sqrt{B}\nabla\cdot\left[\frac{\eta_i}{B}\mathbf{b}\mathbf{b}\cdot\nabla\left(\sqrt{B}V_{||i}\right)\right]

The ion parallel viscosity is

.. math::

   \eta_i = \frac{4}{3} 0.96 p_i \tau_i

The choice of collision frequency is set by the flag `viscosity_collisions_mode`: `multispecies` uses
all available collision frequencies involving the chosen species, while `braginskii` uses only
ii collisions. The default is `multispecies` and it is recommended when solving
more than one ion. If you are solving for a single ion and want to recover Braginskii,
use the `braginskii` mode.


If the `perpendicular` option is set:

.. code-block:: ini

   [braginskii_ion_viscosity]
   perpendicular = true # Include perpendicular flows

Then the ion scalar viscous pressure is calculated as:

.. math::

   \Pi_{ci} = \Pi_{ci||} + \Pi_{ci\perp}

where :math:`\Pi_{ci||}` corresponds to the parallel diffusion of momentum above.

.. math::

   \Pi_{ci||} = - 0.96 \frac{2p_i\tau_i}{\sqrt{B}} \partial_{||}\left(\sqrt{B} V_{||i}\right)

The perpendicular part is calculated from:

.. math::

   \begin{aligned}\Pi_{ci\perp} =& 0.96 p_i\tau_i \kappa \cdot \left[\mathbf{V}_E + \mathbf{V}_{di} + 1.61\frac{\mathbf{b}\times\nabla T_i}{B} \right] \\
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


.. note::
   Implementation details: The magnitude of :math:`\Pi_{ci\perp}` and :math:`\Pi_{ci||}` are
   individually limited to be less than or equal to the scalar pressure :math:`Pi` (though can have
   opposite sign). The reasoning is that if these off-diagonal terms become large then the model is
   likely breaking down. Occasionally happens in low-density regions.


.. doxygenstruct:: BraginskiiIonViscosity
   :members:

.. _sec-braginskii_thermal_force:

braginskii_thermal_force
~~~~~~~~~~~~~~~~~~~~~~~~

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

The implementation is in the `BraginskiiThermalForce` class:

.. doxygenstruct:: BraginskiiThermalForce
   :members:


Neutral gas models
------------------

In 1D, neutral transport is currently done through the same components as for plasma, i.e. `evolve_density`,
`evolve_momentum` and `evolve_pressure` with the additional, optional `neutral_parallel_diffusion` component.
In 2D, all of this functionality is implemented in one component called `neutral_mixed` which additionally
has cross-field transport. This discrepancy is due to historical reasons and will be refactored.


.. _sec-neutral_parallel_diffusion:

1D: neutral_parallel_diffusion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This adds diffusion to **all** neutral species (those with no or zero charge),
because it needs to be calculated after the collision frequencies are known.

.. code-block:: ini

   [hermes]
   components = ... , braginskii_collisions, neutral_parallel_diffusion

   [neutral_parallel_diffusion]
   dneut = 1         # Diffusion multiplication factor
   diagnose = true   # This enables diagnostic output for each species


It is intended mainly for 1D simulations, to provide effective parallel
diffusion of particles, momentum and energy due to the projection of
cross-field diffusion:

.. math::

   \begin{aligned}
   \frac{\partial n_n}{\partial t} =& \ldots + \nabla\cdot\left(\mathbf{b}D_n n_n \frac{1}{p_n}{\partial_{||}p_n}\right) \\
   \frac{\partial p_n}{\partial t} =& \ldots + \nabla\cdot\left(\mathbf{b}D_n p_n \frac{1}{p_n}\partial_{||}p_n\right) + \frac{2}{3}\nabla\cdot\left(\mathbf{b}\kappa_n \partial_{||}T_n\right) \\
   \frac{\partial}{\partial t}\left(m_nn_nv_{||n}\right) =& \ldots + \nabla\cdot\left(\mathbf{b}D_n m_n n_nv_{||n} \frac{1}{p_n} \partial_{||}p_n\right) + \nabla\cdot\left(\mathbf{b}\eta_n \partial_{||}v_{||n}\right)
   \end{aligned}

The diffusion coefficient is in :math:`m^2/s` and is calculated as

.. math::

   D_n = \left(\frac{B}{B_{pol}}\right)^2 \frac{eT_n}{m_{n} \nu}

where ``m_{n}`` is the neutral species mass in kg and :math:`\nu` is the collision
frequency (by default, this sums up all of the enabled neutral collisions from
the collisions component as well as the charge exchange rate).
The factor :math:`B / B_{pol}` is the projection of the cross-field
direction on the parallel transport, and is the ``dneut`` input setting. Currently, the recommended
use case for this component is to represent the neutrals diffusing orthogonal to the target wall, and
it is recommended to set ``dneut`` according to the field line pitch at the target.

.. doxygenstruct:: NeutralParallelDiffusion
   :members:

.. _sec-neutral_mixed:

2D/3D: neutral_mixed
~~~~~~~~~~~~~~~~~~~~~~~~~~

Equations
^^^^^^^^^

The below describes the `neutral_mixed` component used for 2D and 3D simulations. See :ref:`sec-neutral_mixed_bc` for the
corresponding boundary conditions and :ref:`sec-neutral_boundary` on enabling wall energy losses due to neutral reflection.

The `neutral_mixed` component solves fluid equations along :math:`y`
(parallel to the magnetic field), and uses diffusive transport in :math:`x`
and :math:`z`.  It was adopted from the approach used in UEDGE and this [M.V. Umansky, J.N.M (2003)]. The Hermes-3 approach
is more advanced in having a separate neutral pressure equation, similar to the
new AFN (Advanced Fluid Neutral) model in SOLPS-ITER [N. Horsten, N.F. (2017)].

.. math::

   \begin{aligned}

   \frac{\partial n_n}{\partial t} =& -\nabla\cdot\left(n_n\mathbf{b}v_{\parallel, n} + n_n\mathbf{v}_{\perp n}\right) \\
         &    + S \\
   \frac{\partial}{\partial t}\left(m_nn_nv_{\parallel, n}\right) =& -m_n \nabla\cdot\left(n_n v_{\parallel, n} (\mathbf{b}v_{||n} + \mathbf{v}_{\perp n}\right)) \\
         &    - \nabla_{\parallel}p_n \\
         &    + \nabla \cdot (\eta_{n,\perp} \nabla_{\perp} v_{\parallel n} + \mathbf{b} \eta_{n,\parallel} \nabla_{\parallel} v_{\parallel n} ) \\
         &    + F \\
   \frac{\partial p_n}{\partial t} =& -\frac{5}{3} \nabla\cdot\left(p_n\mathbf{b}v_{\parallel, n} +  p_n\mathbf{v}_{\perp n}\right) \\
         &    + \frac{2}{3}v_{\parallel, n} \nabla_{\parallel} p_n \\
         &    + \frac{2}{3} \nabla\cdot\left(\kappa_{n,\perp} \nabla_\perp T_n + \mathbf{b} \kappa_{n,\parallel} \nabla_{\parallel} T_n\right) \\
         &    - \frac{2}{3} v_{\parallel,n} \nabla \cdot \left( \eta_{n,\perp} \nabla_{\perp} v_{\parallel n} + \mathbf{b} \eta_{n,\parallel} \nabla_{\parallel} v_{\parallel n} \right) \\
         &    + \frac{2}{3}E \\

   \end{aligned}

Where for the density equation, the first row of terms contains the parallel and perpendicular
advection and the second row the particle sources. In the parallel momentum equation, the first row of terms
features parallel and perpendicular advection of parallel momentum. This is followed by the compression term
and the perpendicular and parallel viscosity (diffusion of parallel momentum) as well as the momentum source term.
In the pressure equation, the first row contains the parallel and perpendicular advection of pressure. This is followed
by the compression term, the perpendicular and parallel conduction (diffusion of temperature) and perpendicular and parallel
viscous heating, finally followed by the energy sources.

While parallel momentum is evolved and is exchanged with the plasma parallel momentum, the advection of momentum is neglected in the perpendicular direction.
Instead, perpendicular advection is achieved through a pressure diffusion approach, where the pressure gradient is balanced by frictional forces and inertial terms are dropped.
This is similar to Fickian diffusion with the pressure gradient replacing the density gradient as the flow driver. This approach is similar to that taken in nuclear fission neutronic
transport modelling and is typical of fluid neutral models in other codes.

The perpendicular neutral particle flux is:

.. math::
   \begin{aligned}
   \Gamma_{\perp,n_n} =& n_n \mathbf{v}_{\perp},
   \end{aligned}

with the perpendicular velocity :math:`v_{\perp}` calculated as:

.. math::

   \begin{aligned}
   v_{\perp} =& -D_n \frac{1}{p_n} \nabla_{\perp} p_n = -D_n \nabla_{\perp} \ln p_n.
   \end{aligned}

In the code, :math:`\frac{1}{P_n} \nabla_{\perp}P_n` is represented as :math:`\nabla_{\perp} \ln P_n` for numerical reasons.

The unlimited diffusion coefficient is defined as:

.. math::

   \begin{aligned}
   D_{n,unlim} =& \frac{v_{th,n}^{2}}{\nu_{n}}
      = \frac{T_n}{m_n \nu_n} \\
   \end{aligned}

Here :math:`v_{th,n} = \sqrt{T_n / m_n}` is the neutral thermal speed used in
the diffusion coefficient and :math:`\nu_{n,eff}` is the neutral
collision frequency after flooring. A soft floor is applied to prevent
:math:`D_{n,unlim}` from becoming unphysically large at
low collisionality, ensuring that :math:`\nu_n` never goes below a pseudo collision frequency
calculated from the maximum neutral mean free path :math:`l_{\max}`:

.. math::
   \begin{aligned}
   \nu_{mfp} =& \frac{v_{th,n}}{l_{\max}} \\
   \nu_{n} =& \nu_{n}
      + \nu_{mfp}\exp\left(-\frac{\nu_{n}}{\nu_{mfp}}\right) .
   \end{aligned}

Here :math:`l_{\max}` is set by the ``neutral_lmax`` parameter with a default of 1 m.
This is primarily a regularisation: it prevents division by
zero when collision processes are disabled or very weak. In the limit
:math:`\nu_{n} \ll \nu_{mfp}`, the effective collision
frequency tends to :math:`\nu_{mfp}`; in the opposite limit it tends
to the physical collision frequency.

The value of :math:`\nu_{n}` depends on the selected collision mode. When the
``afn`` mode is selected using the
``diffusion_collisions_mode`` setting, the collision frequency is built from
charge exchange, ionisation, and neutral-neutral collisions. When
``multispecies`` mode is selected, Hermes-3 uses charge exchange and all collisions
enabled for neutrals (set in the `braginskii_collisions` component), but ionisation
is neglected.  ``afn`` is the recommended mode, and is based on the Advanced Fluid Neutral
model in SOLPS-ITER.

Flux limiting
^^^^^^^^^^^^^

The diffusive closure places no bound on the transport speed. As the neutral
mean free path grows, the collision frequency falls and the transport
coefficients tend towards infinity. This unphysical behaviour is mitigated
by flux limitation where a flux cap is calculated from a user-specified fraction
of the free-streaming flux.
Hermes-3 limits all three diffusive transport channels: perpendicular
particle diffusion :math:`D_n`, heat conduction :math:`\kappa_n` and viscosity
:math:`\eta_n`.

Each channel is treated in the same four steps:

1. The unlimited coefficient is calculated from the collisionality, as described
   above.
2. The maximum coefficient according to the desired free-streaming fraction
   is calculated.
3. Optionally, an additional, explicit user-set diffusion coefficient cap
   is blended with the free-streaming maximum. The model can also run with a
   purely explicit cap and no free-streaming based limit.
4. The unlimited coefficient and the maximum are blended into the final
   coefficient using a smooth limiter function.

Each flux is the product of a coefficient and a gradient, so step 2 is a
division by the local gradient. That gradient is first regularised using the
function :math:`R`, which is described in the next section:

.. math::

   \begin{aligned}
   D_{n,\max} =& \alpha_{D} \frac{1}{4} \bar{v}_{n}
       \;\Big/\; R\!\left(\nabla_{\perp}\ln p_n\right) \\
   \kappa_{n,\max}^{\perp,\parallel} =& \alpha_{\kappa}^{\perp,\parallel}
       \frac{1}{2} \bar{v}_{n} n_n
       \;\Big/\; R\!\left(\nabla_{\perp,\parallel} T_n / T_n\right) \\
   \eta_{n,\max}^{\perp,\parallel} =& \alpha_{\eta}^{\perp,\parallel} p_n
       \;\Big/\; R_{\eta}\!\left(\nabla_{\perp,\parallel} v_{\parallel n}\right) \\
   \end{aligned}

Here :math:`\bar{v}_{n} = \sqrt{8T_n / \pi m_n}` is the mean speed of a
non-drifting Maxwellian [Stangeby eq. 2.21, p.67]. The numerators are the
free-streaming fluxes of each quantity: :math:`\frac{1}{4} n_n \bar{v}_n` for
particles [Stangeby, under eq. 2.24, p.67], :math:`\frac{1}{2} p_n \bar{v}_n`
for heat, and :math:`p_n` for parallel momentum. The density factor cancels in
:math:`D_{n,\max}` because the perpendicular particle flux is
:math:`\Gamma_{\perp} = -D_n n_n \nabla_{\perp} \ln p_n`.

Conduction and viscosity are limited separately in the parallel and
perpendicular directions, each using the corresponding component of its
gradient. Only the perpendicular direction is
limited for particle transport, because parallel particle transport is
advective rather than diffusive.

Regularising the gradient
^^^^^^^^^^^^^^^^^^^^^^^^^

The gradient in the denominator is smoothly clamped to the range
:math:`[g_{\min}, g_{\max}]` before the division:

.. math::

   \begin{aligned}
   R(g) =& \sqrt{\frac{g^{2} g_{\max}^{2}}{g^{2} + g_{\max}^{2}} + g_{\min}^{2}}
   \end{aligned}

This has three limits:

- :math:`|g| \ll g_{\min}`: :math:`R \rightarrow g_{\min}`. The floor prevents
  division by zero where the gradient vanishes. It is applied as a sum of
  squares so that :math:`R` remains smooth at :math:`g = 0` instead of
  developing a kink. The cost is a small additional reduction of the flux limit
  in shallow gradients.
- :math:`g_{\min} \ll |g| \ll g_{\max}`: :math:`R \rightarrow |g|`, so the
  gradient passes through unmodified.
- :math:`|g| \gg g_{\max}`: :math:`R \rightarrow g_{\max}`. The ceiling stops
  the maximum coefficient from becoming arbitrarily small at very steep
  gradients, which improves robustness where the limiter is saturated. In the
  saturated limit the neutral transport depends on the free-streaming speed
  without a dependence on a gradient, making it behave like an advective term.
  This makes the central difference scheme implemented in its operator unstable
  due to the risk of checkerboarding. The gradient ceiling ensures that the term
  can never become fully advective and improves simulation robustness during transients.
  The default value is set so that it has very little impact in steady state.

Two sets of floor and ceiling parameters are used, because the gradients do not
all have the same units. The particle and conduction limiters both divide by an
inverse gradient length in :math:`m^{-1}` and share one pair, denoted :math:`R`
above. The viscosity limiter divides by a velocity gradient in :math:`s^{-1}`
and has its own pair, denoted :math:`R_{\eta}`. All four parameters are set in
SI units.

.. list-table::
   :header-rows: 1
   :widths: 26 14 12 24 24

   * - Gradient
     - Used in
     - Units
     - Floor option (default)
     - Ceiling option (default)
   * - :math:`\nabla_{\perp}\ln p_n`
     - :math:`D_{n,\max}`
     - :math:`m^{-1}`
     - ``limiter_gradient_floor`` (10)
     - ``limiter_gradient_ceiling`` (100)
   * - :math:`\nabla_{\perp,\parallel} T_n / T_n`
     - :math:`\kappa_{n,\max}`
     - :math:`m^{-1}`
     - ``limiter_gradient_floor`` (10)
     - ``limiter_gradient_ceiling`` (100)
   * - :math:`\nabla_{\perp,\parallel} v_{\parallel n}`
     - :math:`\eta_{n,\max}`
     - :math:`s^{-1}`
     - ``limiter_gradient_floor_eta`` (9.5788e4)
     - ``limiter_gradient_ceiling_eta`` (1e12)

Raising a floor or lowering a ceiling improves robustness but changes the
answer, so both should be checked for convergence.

Applying the limit
^^^^^^^^^^^^^^^^^^

The final coefficient is a smooth blend of the unlimited coefficient and the
maximum. Writing this for :math:`D_n`, and with the same form used for
:math:`\kappa_n` and :math:`\eta_n`:

.. math::

   \begin{aligned}
   D_n = D_{n,unlim}
      \left[1 + \left(\frac{D_{n,unlim}}{D_{n,\max}}\right)^\gamma
      \right]^{-1/\gamma},
   \end{aligned}

where :math:`\gamma` is set by ``flux_limiter_sharpness``. The default
:math:`\gamma = 1` gives the harmonic mean

.. math::

   \begin{aligned}
   D_n = \frac{D_{n,unlim}D_{n,\max}}
              {D_{n,unlim} + D_{n,\max}} ,
   \end{aligned}

and larger values of :math:`\gamma` make the transition sharper.

The flux limit fraction :math:`\alpha` also acts as a switch on its channel:

- :math:`\alpha < 0` disables the free-streaming limiter, so the unlimited
  coefficient is used unless an explicit cap is set.
- :math:`\alpha = 0` sets the coefficient to zero, switching the channel off.
  Note that with ``combined_limiters = true``, setting ``flux_limit = 0``
  switches off diffusion, conduction and viscosity together.

Explicit coefficient caps
^^^^^^^^^^^^^^^^^^^^^^^^^

Each channel can also be capped by an explicit, user-set maximum coefficient,
using ``diffusion_limit``, ``conduction_limit`` and ``viscosity_limit``. These
are set in SI units and are negative by default, which disables them.

If both the flux limit fraction and the explicit cap are set for a channel, the
maximum coefficient is the harmonic mean of the two. If only the explicit cap is
set, the maximum coefficient is that value. If neither is set, the channel is
not limited at all.

Combined and separate limiters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The conduction and viscosity coefficients are related to the diffusion
coefficient by the Chapman-Enskog approximation:

.. math::

   \begin{aligned}
   \kappa_{n} =& \frac{5}{2} D_n n_n \\
   \eta_{n} =& \frac{2}{5} m_n \kappa_{n} \\
   \end{aligned}

This relation always holds for the unlimited coefficients. Whether it also holds
for the limited coefficients is controlled by ``combined_limiters``.

With ``combined_limiters = true``, which is the default, only :math:`D_n` is
limited, and the limited :math:`\kappa_n` and :math:`\eta_n` are calculated from
the limited :math:`D_n` using the relation above. The limiting is then isotropic
and is driven by the pressure gradient alone. This is cheaper and more robust,
but less physically accurate. In this mode, the conduction and viscosity limiter
options have no effect, and setting any of ``flux_limit_cond_perp``,
``flux_limit_cond_par``, ``flux_limit_visc_perp``, ``flux_limit_visc_par``,
``conduction_limit`` or ``viscosity_limit`` raises an exception rather than
being silently ignored.

With ``combined_limiters = false``, each channel is limited by its own
free-streaming flux and its own gradient, separately in the parallel and
perpendicular directions, as described above.

Flux limiter options
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 18 52

   * - Option
     - Default
     - Meaning
   * - ``combined_limiters``
     - ``true``
     - Derive the :math:`\kappa_n` and :math:`\eta_n` limits from :math:`D_n`
   * - ``flux_limit``
     - 0.5
     - :math:`\alpha_{D}`, fraction of free-streaming particle flux
   * - ``flux_limit_cond_perp``
     - ``flux_limit``
     - :math:`\alpha_{\kappa}^{\perp}`, fraction of free-streaming heat flux
   * - ``flux_limit_cond_par``
     - ``flux_limit_cond_perp``
     - :math:`\alpha_{\kappa}^{\parallel}`
   * - ``flux_limit_visc_perp``
     - ``flux_limit``
     - :math:`\alpha_{\eta}^{\perp}`, fraction of free-streaming momentum flux
   * - ``flux_limit_visc_par``
     - ``flux_limit_visc_perp``
     - :math:`\alpha_{\eta}^{\parallel}`
   * - ``flux_limiter_sharpness``
     - 1.0
     - :math:`\gamma`, sharpness of the blend between unlimited and maximum
   * - ``diffusion_limit``
     - -1
     - Explicit cap on :math:`D_n` in :math:`m^{2} s^{-1}`
   * - ``conduction_limit``
     - -1
     - Explicit cap on :math:`\kappa_n` in :math:`m^{-1} s^{-1}`
   * - ``viscosity_limit``
     - -1
     - Explicit cap on :math:`\eta_n` in :math:`Pa\,s`
   * - ``limiter_gradient_floor``
     - 10
     - :math:`g_{\min}` for :math:`D_n` and :math:`\kappa_n`, in :math:`m^{-1}`
   * - ``limiter_gradient_ceiling``
     - 100
     - :math:`g_{\max}` for :math:`D_n` and :math:`\kappa_n`, in :math:`m^{-1}`
   * - ``limiter_gradient_floor_eta``
     - 9.5788e4
     - :math:`g_{\min}` for :math:`\eta_n`, in :math:`s^{-1}`
   * - ``limiter_gradient_ceiling_eta``
     - 1e12
     - :math:`g_{\max}` for :math:`\eta_n`, in :math:`s^{-1}`

Setting a flux limit fraction or an explicit cap below zero disables it.

Flux limiter diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^

Setting ``diagnose = true`` writes the coefficients at each stage of the
limiting to the output. The names below are for a species ``d``.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Diagnostic
     - Meaning
   * - ``Dnnd_unlimited``, ``kappa_d_unlimited``, ``eta_d_unlimited``
     - Coefficients from the collisionality, before any limiting
   * - ``Dnnd_max``, ``kappa_d_max_perp``, ``kappa_d_max_par``,
       ``eta_d_max_perp``, ``eta_d_max_par``
     - Maximum coefficients, after any explicit cap is blended in
   * - ``Dnnd``, ``kappa_d_perp``, ``kappa_d_par``, ``eta_d_perp``,
       ``eta_d_par``
     - Final coefficients used in the equations

The maximum coefficients are only written if the corresponding limiter is
enabled. With ``combined_limiters = true``, the conduction and viscosity
maximums are derived from ``Dnnd_max`` and are written for diagnostic purposes
only.


.. doxygenstruct:: NeutralMixed
   :members:

2D/3D: neutral_full_velocity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This model evolves the equations for a neutral fluid, assuming
axisymmetry (constant in :math:`Z`), for the density :math:`n_n`,
velocity :math:`\mathbf{v}_n` and pressure :math:`p_n`.

.. math::

   \begin{aligned}
   \frac{\partial n_n}{\partial t} =& -\nabla\cdot\left(n_n\mathbf{v}_n\right) \nonumber \\
   \frac{\partial \mathbf{v}_n}{\partial t} =& - \mathbf{v}_n\cdot\nabla\mathbf{v}_n -\frac{1}{n_n m_n}\nabla p_n + \nabla\cdot\left(\nu \nabla\mathbf{v}\right)\\
   \frac{\partial p_n}{\partial t} =& -\gamma \nabla\cdot\left(p_n\mathbf{v}_n\right) + \left(\gamma - 1\right)\mathbf{v}_n\cdot\nabla p_n + \nabla\cdot\left(n_n \chi_n \nabla T_n\right) \nonumber
   \end{aligned}

where the adiabatic index :math:`\gamma` and dissipation parameters
:math:`\nu` (kinematic viscosity) and :math:`\chi` (thermal
conduction) are constants set in the options:

.. code-block:: ini

   [d]
   type = neutral_full_velocity

   adiabatic_index = 5./3 # Ratio of specific heats
   viscosity = 1.0   # Kinematic viscosity [m^2/s]
   conduction = 1.0  # Heat conduction [m^2/s]

The contravariant components of :math:`\mathbf{v}_n` are evolved in
the same :math:`\left(x,y,z\right)` field-aligned coordinate system as
the plasma.  To evaluate the nonlinear advection term, whilst avoiding
the use of noisy Christoffel symbols coming from derivatives of basis
vectors, these components are transformed into
:math:`\left(R,Z,\phi\right)` cylindrical coordinates, advected, then
transformed back. This is done using matrices which are calculated in
the initialisation stage by finite differences of the input mesh:

.. math::

   \begin{aligned}
   \left(\begin{array}{c}
   \nabla R \\
   \nabla Z\end{array}\right) =& \left(\begin{array}{cc}
   \frac{\partial R}{\partial x} & \frac{\partial R}{\partial y} \\
   \frac{\partial Z}{\partial x} & \frac{\partial Z}{\partial y}\end{array}\right)\left(\begin{array}{c}
   \nabla x \\
   \nabla y\end{array}\right) \\
   =& \left(\begin{array}{cc}
   \texttt{Urx} & \texttt{Ury} \\
   \texttt{Uzx} & \texttt{Uzy} \end{array}\right)\left(\begin{array}{c}
   \nabla x \\
   \nabla y\end{array}\right)
   \end{aligned}

These components are calculated by finite differences of the `Rxy` and
`Zxy` arrays in the input, then adjusted to match the given values of
`hthe` and `Bpxy`:

.. math::

   \sqrt{\left(\frac{\partial R}{\partial y}\right)^2 + \left(\frac{\partial R}{\partial y}\right)^2} = h_\theta


.. math::

   \sqrt{\left(\frac{\partial R}{\partial x}\right)^2 + \left(\frac{\partial R}{\partial x}\right)^2} = 1 / \left(R B_\theta\right)


(Note that this second equality only works if :math:`x` and :math:`y`
are orthogonal).

This matrix is then inverted, to give:

.. math::

   \begin{aligned}
   \left(\begin{array}{c}
   \nabla x \\
   \nabla y\end{array}\right) =& \left(\begin{array}{cc}
   \texttt{Txr} & \texttt{Tyr} \\
   \texttt{Txz} & \texttt{Tyz} \end{array}\right)\left(\begin{array}{c}
   \nabla R \\
   \nabla Z\end{array}\right)
   \end{aligned}

The components of :math:`\mathbf{v}_n` are evolved in contravariant form:

.. math::

   \mathbf{v}_n = v^x \mathbf{e}_x + v^y \mathbf{e}_y + v^z \mathbf{e}_z

These components are stored in the output. In the RHS function the
velocity is converted to covariant form:

.. math::

   \mathbf{v}_n = v_x \nabla x + v_y \nabla y + v_z \nabla z

which is then transformed to :math:`v_r`, :math:`v_Z` and :math:`v_\phi`:

.. math::

   \begin{aligned}
   v_r =& \mathbf{v}_n \cdot \nabla R = \frac{\partial x}{\partial R} v_x + \frac{\partial y}{\partial R} v_y \\
   v_Z =& \mathbf{v}_n \cdot \nabla Z = \frac{\partial x}{\partial Z} v_x +  \frac{\partial y}{\partial Z} v_y \\
   v_\phi =& \mathbf{v}_n \cdot \hat{\phi} = v_z / \left(\sigma_{Bpol} R\right)
   \end{aligned}

which are implemented as

.. code-block:: c++

   Field2D vr = Txr * Vn2D.x + Tyr * Vn2D.y;
   Field2D vz = Txz * Vn2D.x + Tyz * Vn2D.y;
   Field2D vphi = Vn2D.z / (sigma_Bp * Rxy);

These components are then advected as scalars for the
:math:`\mathbf{v}_n\cdot\nabla\mathbf{v}_n` term, and are diffused for
the :math:`\nabla\cdot\left(\mu \nabla\mathbf{v}\right)` kinematic
viscosity.

The advection of momentum :math:`\mathbf{v}\cdot\nabla\mathbf{v}` is
controlled by these settings:

#. `momentum_advection` is ``false`` by default, disabling this
   nonlinear advection term. This keeps the inertia in the time
   derivative, but neglects the neutral dynamic pressure in the
   momentum balance.

#. `toroidal_flow` is ``true`` by default, which includes the toroidal
   (:math:`z`) component of the neutral flow. Importantly, this allows
   the parallel and poloidal flows to evolve independently: The
   parallel flow can follow the plasma towards the target, while the
   poloidal flow can be away from the target.

#. `curved_torus` is ``true`` by default, and is only active when both
   `momentum_advection` and `toroidal_flow` are enabled. Neutrals
   travel in straight lines in real space, so toroidal flow is
   converted to radial flow. This appears in the :math:`v_r` and
   :math:`v_\phi` equations due to a combination of the radial
   centrifugal force and conservation of toroidal angular momentum.

Flow perpendicular to the magnetic field is damped by collisions
e.g. CX reactions with the plasma. The steady-state flow is therefore
a balance between the pressure gradient (including dynamic pressure if `momentum_advection` is enabled),
and this friction. The neutral velocity perpendicular to the magnetic field is:

.. math::

   \begin{aligned}
   \mathbf{v}_{n\perp} =& \mathbf{v}_{n} - \mathbf{b}\mathbf{b}\cdot\mathbf{v}_{n} \\
   =& \mathbf{v}_{n} - \mathbf{e}_y\frac{v_{ny}}{g_{yy}} \\
   =& \mathbf{v}_{n} - \left(\nabla y + \frac{g_{yz}}{g_{yy}}\nabla z\right)v_{ny} \\
   =& v_{nx}\nabla x + \left(v_{nz} - \frac{g_{yz}}{g_{yy}}v_{ny}\right)\nabla z
   \end{aligned}

At boundaries neutral thermal energy is lost at a rate controlled by
the option

.. code-block:: ini

   neutral_gamma = 5./4

This sets the flux of power to the wall to:

.. math::

   q = \gamma n_n T_n c_s

Currently this is only done at target boundaries, not radial
boundaries.

Drifts and transport
--------------------

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

Two forms are available, which are implemented differently for density, momentum, and pressure equations. In the density equation, form 0 uses the diamagnetic velocity perpendicular to b and the gradient of P;
at the boundaries this velocity is perpendicular to the boundary. Form 1 uses the magnetic gyro-center drifts, which are mostly vertical;
at the boundaries this form produces a flow through the boundary.
Forms 0 and 1 are analytically equivalent and should give the same result away from boundaries,
but form 0 doesn't produce flows through boundaries. This is an approach that UEDGE uses to avoid unphysical boundary flows.


However, Form 1 is nice because the flow velocity depends on the temperature, not the pressure gradient.
This usually makes it better behaved numerically. To make the most of both, the `diamagnetic_drift` component allows the forms to be mixed
using the ``diamag_form`` setting. For example, the :code:`tcv-x21` example blends it such that form 0 is at the boundary:


.. code-block:: ini

   [diamagnetic_drift]
   diamag_form = x * (1 - x)  # 0 = gradient; 1 = divergence


A table of the two forms used in Hermes-3, and the corresponding terms in `Simakov & Catto <https://doi.org/10.1063/1.1623492>`_ is shown below, where :math:`\mathbf{C}=\nabla\times\left(\frac{\mathbf{b}}{B}\right)` is the curvature vector. Instead of the diamagnetic velocity, the whole terms associated are shown. The difference among the forms is the divergence of a curl, which vanishes. The diamagnetic velocity :math:`\mathbf{v}_{dia}` is defined above. Notice that Simakov & Catto used Gaussian units, but Hermes-3 uses SI units.

.. list-table::
   :header-rows: 1
   :widths: 10 20 20 35

   * -
     - Form 0
     - Form 1
     - Simakov & Catto
   * - Density
     - :math:`\mathbf{C} \cdot \nabla\left(\dfrac{p}{q}\right)`
     - :math:`\nabla \cdot (n \mathbf{v}_{dia})`
     - Eq. (51): :math:`\dfrac{c}{q}\left(\nabla \times \dfrac{\mathbf{b}}{B}\right) \cdot \nabla p`
   * - Momentum
     - :math:`\mathbf{C} \cdot \nabla\left(\dfrac{mnv_\parallel T}{q}\right)`
     - :math:`\nabla\cdot (mnv_\parallel \mathbf{v}_{dia})`
     - Eq. (64): :math:`\nabla\cdot \left(\dfrac{1}{\Omega} \mathbf{b} \times \nabla(p v_\parallel)\right)`
   * - Pressure
     - :math:`\dfrac{5}{2}\mathbf{C} \cdot \nabla\left(\dfrac{pT}{q}\right)`
     - :math:`\dfrac{5}{2}\nabla\cdot (p \mathbf{v}_{dia})`
     - Eq. (56): :math:`\nabla\cdot\left(\dfrac{5}{2 m \Omega} \mathbf{b} \times \nabla(p T)\right)`

\* Eq.(64) in Simakov & Catto is derived for ion parallel momentum, but it is also applicable to electrons since it comes from the gyro-viscosity and the mass factors of :math:`m_i` or :math:`m_e` cancel out.


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

See the :file:`examples/blob2d-vpol` example, which contains:

.. code-block:: ini

   [hermes]
   components = e, vorticity, sheath_closure, polarisation_drift

   [polarisation_drift]
   diagnose = true

Setting ``diagnose = true`` saves `DivJ` to the dump files with the divergence of all
currents except polarisation, and `phi_pol` which is the polarisation flow potential.

.. doxygenstruct:: PolarisationDrift
   :members:

Stellarator cross-field transport: binormal_stpm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This adds a term to **all** species which includes the effects of cross-field
drifts following the stellarator two point model:
`Y. Feng et al., Plasma Phys. Control. Fusion 53 (2011) 024009 <http://dx.doi.org/10.1088/0741-3335/53/2/024009>`_

.. code-block:: ini

   [hermes]
   components = ... , binormal_stpm

   [binormal_stpm]
   D = 1         # [m^2/s]  Density diffusion coefficient
   chi = 3       # [m^2/s]  Thermal diffusion coefficient
   nu = 1        # [m^2/s]  Momentum diffusion coefficient

   Theta = 1e-3  # Field line pitch

It is intended only for 1D simulations, to provide effective parallel
diffusion of particles, momentum and energy due to the projection of
cross-field diffusion:

.. math::

   \begin{aligned}
   \frac{\partial N}{\partial t} =& \ldots + \nabla\cdot\left(\mathbf{b}\frac{D}{\Theta}\partial_{||}N\right) \\
   \frac{\partial P}{\partial t} =& \ldots + \frac{2}{3}\nabla\cdot\left(\mathbf{b}\frac{\chi}{\Theta} N\partial_{||}T\right) \\
   \frac{\partial}{\partial t}\left(NV\right) =& \ldots + \nabla\cdot\left(\mathbf{b}\frac{\nu}{\Theta} \partial_{||}NV\right)
   \end{aligned}

The diffusion coefficients `D`, `\chi` and `\nu` and field line pitch `\Theta` are prescribed in the input file.


.. doxygenstruct:: BinormalSTPM
   :members:


Tokamak cross-field transport: anomalous_diffusion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adds cross-field diffusion of particles, momentum and energy to a species.

.. code-block:: ini

   [hermes]
   components = e, ...

   [e]
   type = evolve_density, evolve_momentum, evolve_pressure, anomalous_diffusion

   anomalous_D = 1.0   # Density diffusion [m^2/s]
   anomalous_chi = 0,5 # Thermal diffusion [m^2/s]
   anomalous_nu = 0.5  # Kinematic viscosity [m^2/s]

Anomalous diffusion coefficients can be functions of `x` and `y`.  The
coefficients can also be read from the mesh input file: If the mesh
file contains `D_` + the name of the species, for example `D_e` then
it will be read and used as the density diffusion coefficient.
Similarly, `chi_e` is the thermal conduction coefficient, and `nu_e`
is the kinematic viscosity. All quantities should be in SI units of
m^2/s.  Values that are set in the options (as above) override those
in the mesh file.

The sources of particles :math:`S`, momentum :math:`F` and energy
:math:`E` are calculated from species density :math:`N`, parallel
velocity :math:`V` and temperature :math:`T` using diffusion
coefficients :math:`D`, :math:`\chi` and :math:`\nu` as follows:

.. math::

   \begin{aligned}
   S =& \nabla\cdot\left(D \nabla_\perp N\right) \\
   F =& \nabla\cdot\left(m V D \nabla_\perp N\right) + \nabla\cdot\left(m N \nu \nabla_\perp V\right)\\
   E =& \nabla\cdot\left(\frac{3}{2}T D \nabla_\perp N\right) + \nabla\cdot\left(N \chi \nabla_\perp T\right)
   \end{aligned}

Note that particle diffusion is treated as a density gradient-driven flow
with velocity :math:`v_D = -D \nabla_\perp N / N`.

.. doxygenstruct:: AnomalousDiffusion
   :members:


Electromagnetic fields
----------------------

These are components which calculate the electric and/or magnetic
fields.

.. _sec-vorticity:

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


Kinematic viscosity can be included by setting e.g. ``viscosity = 0.1`` in SI units (m^2/s).
This adds a diffusion of vorticity and corresponding ion heating.
The viscous friction force in this simplified operator is

.. math::

   \mathbf{F}_\nu = - \nu B \mathbf{b} \times \nabla \Omega

This gives rise to a drift and current with divergence:

.. math::

   \begin{aligned}\nabla\cdot\mathbf{J_{\nu}} =& \nabla\cdot\left[\frac{\mathbf{b}\times\mathbf{F}_\nu}{B}\right] \\
   =& \nabla\cdot\left[\nu \nabla_\perp \Omega\right]\end{aligned}

Viscous heating is calculated using the work done by a fluid velocity
consistent with the Boussinesq approximation:

.. math::

   \mathbf{u} = \frac{\mathbf{b}\times\nabla\Phi}{B}

where the generalized potential is

.. math::

   \Phi = \phi + \hat{p} / \overline{n}

The work done is

.. math::

   \mathbf{F}_\nu\cdot\mathbf{u} = -\nu \nabla_\perp\Omega \cdot \nabla_\perp\Phi

This heating is distributed between charged species in proportion to their local mass density.
The properties of the work done can be analysed by writing in terms of a vector :math:`\mathbf{g}`:

.. math::

   \mathbf{g} = \frac{\overline{A}\overline{n}}{B^2}\nabla_\perp\Phi

to write:

.. math::

   \begin{aligned}\mathbf{F}_\nu\cdot\mathbf{u} =& -\nu\frac{B^2}{\overline{A}\overline{n}} \nabla_\perp\left(\nabla\cdot\mathbf{g}\right)\cdot\mathbf{g} \\
   =& \nu\frac{B^2}{\overline{A}\overline{n}} \left[\left(\nabla_\perp\mathbf{g} : \nabla_\perp\mathbf{g}\right) - \nabla\cdot\left(\mathbf{g}\cdot\nabla\mathbf{g}\right)\right]\end{aligned}

The last term is not in general positive definite, so this simple form
of viscosity could in some cases lead to cooling.

.. doxygenstruct:: Vorticity
   :members:

relax_potential
~~~~~~~~~~~~~~~

This component evolves a vorticity equation, similar to the ``vorticity`` component.
Rather than inverting an elliptic equation at every timestep, this component evolves
the potential in time as a diffusion equation.

.. doxygenstruct:: RelaxPotential
   :members:

.. _sec-electromagnetic:

electromagnetic
~~~~~~~~~~~~~~~

**Notes**: When using this module,

1. Set ``sound_speed:alfven_wave=true`` so that the shear Alfven wave
   speed is included in the calculation of the fastest parallel wave
   speed for numerical dissipation.
2. For tokamak simulations use Neumann boundary condition on the core
   and Dirichlet on SOL and PF boundaries by setting
   ``electromagnetic:apar_core_neumann=true`` (this is the default).
3. Set the potential core boundary to be constant in Y by setting
   ``vorticity:phi_core_averagey = true``
4. Magnetic flutter terms must be enabled to be active
   (``electromagnetic:magnetic_flutter=true``).  They use an
   ``Apar_flutter`` field, not the ``Apar`` field that is calculated
   from the induction terms.
5. If using ``vorticity:phi_boundary_relax`` to evolve the radial
   boundary of the electrostatic potential, the timescale
   ``phi_boundary_timescale`` should be set much longer than the
   Alfven wave period or unphysical instabilities may grow from the
   boundaries.

This component modifies the definition of momentum of all species, to
include the contribution from the electromagnetic potential
:math:`A_{||}`.

Assumes that "momentum" :math:`p_s` calculated for all species
:math:`s` is

.. math::

   p_s = m_s n_s v_{||s} + Z_s e n_s A_{||}

which arises once the electromagnetic contribution to the force on
each species is included in the momentum equation. This requires
an additional term in the momentum equation:

.. math::

   \frac{\partial p_s}{\partial t} = \cdots + Z_s e A_{||} \frac{\partial n_s}{\partial t}

This is implemented so that the density time-derivative is calculated using the lowest order
terms (parallel flow, ExB drift and a low density numerical diffusion term).

The above equations are normalised so that in dimensionless quantities:

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

The toroidal variation of density :math:`n_s` must be kept in this
equation.  By default the iterative "Naulin" solver is used to do
this: A fast FFT-based method is used in a fixed point iteration,
correcting for the density variation.

Magnetic flutter terms are disabled by default, and can be enabled by setting

.. code-block:: ini

   [electromagnetic]
   magnetic_flutter = true

This writes an ``Apar_flutter`` field to the state, which then enables perturbed
parallel derivative terms in the ``evolve_density``, ``evolve_pressure``, ``evolve_energy`` and
``evolve_momentum`` components. Parallel flow terms are modified, and parallel heat
conduction.

.. math::

   \begin{aligned}\mathbf{b}\cdot\nabla f =& \mathbf{b}_0\cdot\nabla f + \delta\mathbf{b}\cdot\nabla f \\
   =& \mathbf{b}_0\cdot\nabla f + \frac{1}{B}\nabla\times\left(\mathbf{b}A_{||}\right)\cdot\nabla f \\
   \simeq& \mathbf{b}_0\cdot\nabla f + \frac{1}{B_0}\left[A_{||}\nabla\times\mathbf{b} + \left(\nabla A_{||}\right)\times\mathbf{b}_0\right]\cdot\nabla f \\
   \simeq& \mathbf{b}_0\cdot\nabla f + \frac{1}{B_0}\mathbf{b}_0\times \nabla A_{||} \cdot \nabla f\end{aligned}

.. doxygenstruct:: Electromagnetic
   :members:
