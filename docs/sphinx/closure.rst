.. _sec-closure:

Closure
==========

Hermes-3 currently uses the standard Braginskii closure, but with a selection of 
collision frequencies that can be used: the formulas use all solved collisions
by default in order to have some accounting of multiple ion species, which 
are not considered in standard Braginskii. See the next section for how this can be
changed to reproduce the exact Braginskii closure.

The different parts of the closures are currently implemented in different components,
which are described in the :ref:`sec-equations` section. They all use collision frequencies
calculated in the ``Collisions`` component.

**Conduction:** 
Parallel conduction for any species is in the :ref:`evolve_pressure` and 
:ref:`evolve_energy` species-level components.

**Viscosity:**
Parallel and perpendicular viscosity for ions is in the 
:ref:`ion_viscosity` top-level component. The parallel viscosity for 
electrons is in :ref:`electron_viscosity`.

**Thermal force:**
Thermal force is implemented here: :ref:`thermal_force`.

**Frictional and thermal equilibration:**
both are calculated in the top-level
``collisions`` component described in section :ref:`sec-collisions`.

**Neutral diffusion:**
The parallel projection of diffusion from the wall in 1D
is captured in the :ref:`neutral_parallel_diffusion` top-level component, while 
both parallel Braginskii transport and perpendicular pressure-diffusion for 2D/3D 
are captured in the :ref:`neutral_mixed` species-level component. 
   


Collision frequency selection
~~~~~~~~~~

When configured in ``multispecies``
mode (default), all collision frequencies for solved collisions of a particular 
species are counted, e.g. electron conduction would include ``ee``, ``ei`` and ``en``
collisions, and ion conduction would include ``ii``, ``ie``, ``in`` and ``CX``.
In this way, the closure is similar to that used in UEDGE and has
a way to account for multiple ion species.
Note that any of the collisions can be disabled (see next section).

For code comparison purposes, it can be useful to reproduce exact Braginskii 
closure, e.g. featuring only self-collisions for conduction and viscosity.
This is enabled by changing the mode to ``braginskii``. Note that this will be invalid
for simulations with multiple ion species.

The collision frequency choice for conduction is done under the species header, e.g.:

.. code-block:: ini

   [d+]  # Deuterium ions
   type = (evolve_density, evolve_pressure, evolve_momentum,
         noflow_boundary, upstream_density_feedback)

   conduction_collisions_mode = braginskii

The choice for viscosity is done under the ``[ion_viscosity]`` header due to it
being a top level component, e.g.:

.. code-block:: ini

   [ion_viscosity]
   viscosity_collisions_mode = braginskii

In addition to the parallel closure, there is also a collision frequency choice
for neutral diffusion, which is present both in ``neutral_parallel_diffusion`` (1D) 
and ``neutral_mixed`` (2D). In this case, there are again two choices:
``multispecies`` features all enabled collisions, while ``afn`` selects
only ionisation and charge exchange, consistent with the AFN (Advanced Fluid Neutral)
work in `N. Horsten Nucl. Fusion 57 (11) 116043 (2017) <https://doi.org/10.1088/1741-4326/AA8009>`_.

``neutral_mixed`` is a species level component and should be set under the species header:

.. code-block:: ini

   [d]
   type = neutral_mixed
   diffusion_collisions_mode = afn

While ``neutral_parallel_diffusion`` is a top level component and must be set under its own 
header, e.g.:

.. code-block:: ini

   [neutral_parallel_diffusion]
   diffusion_collisions_mode = afn



.. _sec-collisions:

Collisions component
~~~~~~~~~~

Inputs and ouputs
----------

This top-level component calculates the collision frequencies of all collisional processes
in Hermes-3. These frequencies are then used to calculate the closure terms.
By default, the following collisions are enabled:

.. code-block:: ini

   [collisions]
   electron_ion = true
   electron_electron = true
   electron_neutral = false
   ion_ion = true
   ion_neutral = false
   neutral_neutral = true

``electon_neutral`` collisions are disabled as they are are typically
a very minor contributor, while ``ion_neutral`` collisions are disabled as 
they are already accounted for by charge exchange which is enabled by default.

All of the collision frequencies are added to the state in ``species["collision_frequencies"]``.
They are also available as diagnostics, e.g. ``Kd+e_coll`` is the ion-electron collision
frequency. 

Theory
----------


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

The electron-electron collision time definition follows Braginskii (note that Fitzpatrick uses 
a different definition in his `notes <https://farside.ph.utexas.edu/teaching/plasma/Plasma/node41.html>`_,
these are not consistent with Braginskii):

.. math::
   \nu_{ee} = \frac{ln \Lambda e^4 n_e} { 12 \pi^{3/2} \varepsilon_0^2 m_{e}^{1/2} T_{e}^{3/2} } 

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

  *Note*: These are disabled by default. If enabled, care is needed to
  avoid double-counting collisions in atomic reactions e.g charge-exchange
  reactions.
  
  The cross-section for elastic collisions between charged and neutral
  particles can vary significantly. Here for simplicity we just take
  a value of :math:`5\times 10^{-19}m^2` from the NRL formulary.

- Neutral-neutral collisions

  *Note* This is enabled by default.
  
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