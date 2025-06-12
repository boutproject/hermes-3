.. _sec-closure:

Closure
==========

Hermes-3 currently uses a closure similar to that in UEDGE. It uses
the standard Braginskii equations, but with additional collision frequencies
added due to additional species. For example, Braginskii ion conduction 
would enable collisions from any enabled processes and with any species,
e.g. electron-ion collisions, ion-ion collisions with any other species,
as well as charge exchange and ionisation. The ``[collisions]`` component 
detailed in this section controls which collisional processes are enabled.

This method allows the closure to account for all of the species present
in the simulation, but it is a very simple implementation as Braginskii
was derived for a single main ion species. Work on improved closure
is ongoing.

The closures for conduction and viscosity are located in the equation components
``evolve_pressure`` / ``evolve_energy`` and ``evolve_momentum``. Please
see the Equations section for more details.


Collisions
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