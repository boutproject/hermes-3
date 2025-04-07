.. _sec-reactions:

Reactions
===========

The following content gives some background to the reactions implemented in Hermes-3.
Please see the :ref:`postprocessing` section for related diagnostics.


Reaction basics
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







Transfer channels
~~~~~~~~~~~

Reactions typically convert species from one to another, leading to
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

Adjusting reactions
-----------

The reaction rates can be adjusted by a user-specified arbitrary multiplier. This can be useful for 
the analysis of the impact of individual reactions. The multiplier setting must be placed under the 
neutral species corresponding to the reaction, e.g. under `[d]` when adjusting deuterium ionisation, recombination or charge exchange.
The multiplier for the fixed fraction impurity radiation must be placed under the impurity species header, e.g. under `[ar]` for argon.
This functionality is not yet currently implemented for helium or neon reactions.

+-----------------------+------------------+---------------------------------------+
| Setting               | Specified under  |  Reaction                             |
+=======================+==================+=======================================+
| K_iz_multiplier       | Neutral species  | Ionisation rate                       |
+-----------------------+------------------+---------------------------------------+
| R_ex_multiplier       | Neutral species  | Ionisation (excitation) radiation rate|
+-----------------------+------------------+---------------------------------------+
| K_rec_multiplier      | Neutral species  | Recombination rate                    |
+-----------------------+------------------+---------------------------------------+
| R_rec_multiplier      | Neutral species  | Recombination radiation rate          |
+-----------------------+------------------+---------------------------------------+
| K_cx_multiplier       | Neutral species  | Charge exchange rate                  |
+-----------------------+------------------+---------------------------------------+
| R_multiplier          | Impurity species | Fixed frac. impurity radiation rate   |
+-----------------------+------------------+---------------------------------------+

The charge exchange reaction can also be modified so that the momentum transfer channel is disabled. This can be useful when
testing the impact of the full neutral momentum equation equation compared to purely diffusive neutrals. A diffusive only model 
leads to all of the ion momentum being lost during charge exchange due to the lack of a neutral momentum equation.
Enabling neutral momentum introduces a more accurate transport model but also prevents CX momentum from being lost, which
can have a significant impact on the solution and may be difficult to analyse.
Disabling the momentum transfer channel allows you to study the impact of the improved transport only and is set as:

.. code-block:: ini

   [hermes]
   components = ..., c, ...

   [reactions]
   no_neutral_cx_mom_gain = true


Fixed fraction radiation model
---------------

In the fixed fraction radiation model, the impurity density is assumed to be a constant
fraction of the electron density. The impurity only affects the plasma solution through
radiation, which is calculated using a "cooling curve", which is a function of the radiation
in terms of temperature in units of :math:`Wm^{-3}`. More information on cooling curves
can be found in literature, e.g. `A. Kallenbach PPCF 55(12) (2013) <https://doi.org/10.1088/0741-3335/55/12/124041>`_.



To use this component you can just add it to the list of components and then
configure the impurity fraction:

.. code-block:: ini

   [hermes]
   components = ..., c, ...

   [c]
   type = fixed_fraction_carbon
   fraction = 0.05   # 5% of electron density
   diagnose = true   # Saves Rc (R + section name)

This will create a diagnostic variable named ``Rc`` containing the radiation
source in :math:`Wm^{-3}` where ``c`` is the species name.


Several ADAS rates are provided: nitrogen, neon, argon, krypton, xenon and tungsten.
The component names are ``fixed_fraction_carbon``, ``fixed_fraction_nitrogen``, ``fixed_fraction_neon``,
``fixed_fraction_argon``, ``fixed_fraction_krypton``, ``fixed_fraction_xenon`` and ``fixed_fraction_tungsten``.

Each rate is in the form of a 10 coefficient 
log-log polynomial fit of data obtained using the open source tool `radas <https://github.com/cfs-energy/radas>`_, except
xenon and tungsten that use 15 and 20 coefficients respectively.
The :math:`n {\tau}` parameter representing the density and residence time assumed in the radas 
collisional-radiative model has been set to :math:`1\times 10^{20} \times 0.5ms` based on `David Moulton et al 2017 PPCF 59(6) <https://doi.org10.1088/1361-6587/aa6b13>`_.

Each rate has an upper and lower bound beyond which the rate remains constant. 
Please refer to the source code in `fixed_fraction_radiation.hxx` for the coefficients and bounds used for each rate.

In addition to the above rates, there are three simplified cooling curves for Argon: ``fixed_fraction_argon_simplified1``,
``fixed_fraction_argon_simplified2`` and ``fixed_fraction_argon_simplified3``. They progressively reduce the nonlinearity in the 
rate by taking out the curvature from the slopes, taking out the RHS shoulder and taking out the LHS-RHS asymmetry, respectively.
These rates may be useful in investigating the impact of the different kinds of curve nonlinearities on the solution. 

There is also a very simple carbon radiation function ``fixed_fraction_hutchinson_carbon`` which is
in coronal equilibrium, using a simple formula from `I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994) <https://doi.org/10.1088/0029-5515/34/10/I04>`_:

.. math::

   L\left(T_e\right) = 2\times 10^{-31} \frac{\left(T_e/10\right)^3}{1 + \left(T_e / 10\right)^{4.5}}

which has units of :math:`Wm^3` with :math:`T_e` in eV.

NOTE:
   By default, fixed fraction radiation is disabled in the core region. This represents the fact that
   the impurity will likely be coronal on closed field lines and feature reduced radiation. This 
   can prevent unphysical MARFE-like behaviour in deep detachment. This behaviour can be disabled
   by setting ``no_core_radiation=false`` in the impurity options block.

     
Implemented reactions 
-----------

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

+------------------+----------------------------------------------+
| Reaction         | Description                                  |
+==================+==============================================+
| h + e -> h+ + 2e | Hydrogen ionisation (Amjuel H.4 2.1.5)       |
+------------------+----------------------------------------------+
| d + e -> d+ + 2e | Deuterium ionisation (Amjuel H.4 2.1.5)      |
+------------------+----------------------------------------------+
| t + e -> t+ + 2e | Tritium ionisation (Amjuel H.4 2.1.5)        |
+------------------+----------------------------------------------+
| h + h+ -> h+ + h | Hydrogen charge exchange (Amjuel H.3 3.1.8)  |
+------------------+----------------------------------------------+
| d + d+ -> d+ + d | Deuterium charge exchange (Amjuel H.3 3.1.8) |
+------------------+----------------------------------------------+
| t + t+ -> t+ + t | Tritium charge exchange (Amjuel H.3 3.1.8)   |
+------------------+----------------------------------------------+
| h + d+ -> h+ + d | Mixed hydrogen isotope CX (Amjuel H.3 3.1.8) |
+------------------+----------------------------------------------+
| d + h+ -> d+ + h |                                              |
+------------------+----------------------------------------------+
| h + t+ -> h+ + t |                                              |
+------------------+----------------------------------------------+
| t + h+ -> t+ + h |                                              |
+------------------+----------------------------------------------+
| d + t+ -> d+ + t |                                              |
+------------------+----------------------------------------------+
| t + d+ -> t+ + d |                                              |
+------------------+----------------------------------------------+
| h+ + e -> h      | Hydrogen recombination (Amjuel H.4 2.1.8)    |
+------------------+----------------------------------------------+
| d+ + e -> d      | Deuterium recombination (Amjuel H.4 2.1.8)   |
+------------------+----------------------------------------------+
| t+ + e -> t      | Tritium recombination (Amjuel H.4 2.1.8)     |
+------------------+----------------------------------------------+

In addition, the energy loss associated with the ionisation potential energy cost
as well as the photon emission during excitation and de-excitation during multi-step 
ionisation is calculated using the AMJUEL rate H.10 2.1.5. The equivalent rate
for recombination is H.10 2.1.8.

The code to calculate the charge exchange rates is in
`hydrogen_charge_exchange.[ch]xx`. This implements reaction H.3 3.1.8 from
Amjuel (p43), scaled to different isotope masses and finite neutral
particle temperatures by using the effective temperature (Amjuel p43):

.. math::

   T_{eff} = \frac{M}{M_1}T_1 + \frac{M}{M_2}T_2


The effective hydrogenic ionisation rates are calculated using Amjuel
reaction H.4 2.1.5, by D.Reiter, K.Sawada and T.Fujimoto (2016).
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

Lithium
~~~~~~~

These rates are taken from ADAS ('96 and '89)

+-----------------------+---------------------------------------+
| Reaction              | Description                           |
+=======================+=======================================+
| li + e -> li+ + 2e    | Lithium ionisation                    |
+-----------------------+---------------------------------------+
| li+ + e -> li+2 + 2e  |                                       |
+-----------------------+---------------------------------------+
| li+2 + e -> li+3 + 2e |                                       |
+-----------------------+---------------------------------------+
| li+ + e -> li         | Lithium recombination                 |
+-----------------------+---------------------------------------+
| li+2 + e -> li+       |                                       |
+-----------------------+---------------------------------------+
| li+3 + e -> li+2      |                                       |
+-----------------------+---------------------------------------+
| li+ + h -> li + h+    | Charge exchange with hydrogen         |
+-----------------------+---------------------------------------+
| li+2 + h -> li+ + h+  |                                       |
+-----------------------+---------------------------------------+
| li+3 + h -> li+2 + h+ |                                       |
+-----------------------+---------------------------------------+
| li+ + d -> li + d+    | Charge exchange with deuterium        |
+-----------------------+---------------------------------------+
| li+2 + d -> li+ + d+  |                                       |
+-----------------------+---------------------------------------+
| li+3 + d -> li+2 + d+ |                                       |
+-----------------------+---------------------------------------+
| li+ + t -> li + t+    | Charge exchange with tritium          |
+-----------------------+---------------------------------------+
| li+2 + t -> li+ + t+  |                                       |
+-----------------------+---------------------------------------+
| li+3 + t -> li+2 + t+ |                                       |
+-----------------------+---------------------------------------+

The implementation of these rates is in `ADASLithiumIonisation`,
`ADASLithiumRecombination` and `ADASLithiumCX` template classes:

.. doxygenstruct:: ADASLithiumIonisation
   :members:

.. doxygenstruct:: ADASLithiumRecombination
   :members:

.. doxygenstruct:: ADASLithiumCX
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



