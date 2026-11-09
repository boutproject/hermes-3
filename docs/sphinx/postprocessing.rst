.. _sec-postprocessing:

Post-processing
---------------

Tools
~~~~~~~~~~~

There are two workflows for post-processing Hermes-3 results. The legacy way
is to use the `boutdata <https://github.com/boutproject/boutdata>`_ package, which
can be used to extract the raw variables from the dump files individually, which
then have to be unnormalised to be converted into SI units.

The recommended way is to use `xHermes <https://github.com/boutproject/xhermes>`_, 
which is built upon the Xarray BOUT++ post-processing tool `xBOUT 
<https://github.com/boutproject/xBOUT>`_.
The main advantage of xHermes is that it contains several routines for pre-processing
the data, including automatic unnormalisation and geometry calaculations. The disadvantage
is that it may be slower than ``boutdata`` when running in parallel due to the netCDF
I/O overhead. This can be mitigated by squashing the dump files (see section :ref:`sec-execution-squashing`).
Efforts to improve the performance are ongoing.

Please refer to the xHermes `readme <https://github.com/boutproject/xhermes?tab=readme-ov-file#xhermes>`_
and `examples <https://github.com/boutproject/xhermes/tree/main/examples>`_ for details.

Available diagnostic variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It can be problematic to record all of the diagnostic variables in the documentation
as the code is still in active development. Thankfully, all of the diagnostics are
saved with helpful metadata containing the name, units, normalisation, description
and the origin of the variable from within the code.

This is available in xHermes. You can access it by the following example:

.. code-block:: ini

   ds["Nd"].attrs()

The attributes are:

time_dimension
   This indicates that the variable is time evolving. If you remove this line,
   it will be saved only at the first RHS evaluation.

units
   A string showing the units for post-processing. xHermes picks this up.

conversion
   A float representing the normalisation factor. xHermes picks this up to do
   automatic conversion to SI units.

standard_name, long_name
   We aren't consistent on what should be in each, but they are meant to describe
   the variables in post-processing.

species, source
   The relevant species and component that the diagnostic is coming from

Alternatively, the same attributes can be seen in the source code.


Reaction diagnostics convention
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are five categories of diagnostics. ``K`` is a simple reaction rate (always positive).
``S`` is a source of density and is positive if the ion density is increasing.
``E`` is a source of energy due to energy transfer (e.g. charge exchange) and 
is positive if the ion energy is increasing.
``F`` is a momentum transfer rate and is currently negative if the ion energy is increasing.
``R`` is a loss of energy from the system, e.g. due to radiation. It is always negative.

+------------------+---------------------------+-------------------------+
| Variable prefix  |   Units                   | Description             |
+==================+===========================+=========================+
| K                |   :math:`s^{-1}`          | Reaction rate           |
+------------------+---------------------------+-------------------------+
| S                |   :math:`m^{-3}s^{-1}`    | Density transfer source |
+------------------+---------------------------+-------------------------+
| E                |   :math:`Wm^{-3}`         | Energy transfer source  |
+------------------+---------------------------+-------------------------+
| F                |   :math:`kgm^{-2}s^{-2}`  | Momentum transfer rate  |
+------------------+---------------------------+-------------------------+
| R                |   :math:`Wm^{-3}`         | System loss source      |
+------------------+---------------------------+-------------------------+

Charge exchange can be slightly complicated:

For two species ``a`` and ``b``, the channel ``Fab_cx`` is a source of momentum for species ``a`` due to
charge exchange with species ``b``. There are corresponding sinks for
the products of the charge exchange reaction which are not saved.

For example,reaction ``d + t+ -> d+ + t`` will save the following
forces (momentum sources):
- ``Fdt+_cx`` is a source of momentum for deuterium atoms ``d`` and sink of momentum for deuterium ions ``d+``.
- ``Ft+d_cx`` is a source of momentum for tritium ions ``t+`` and sink of momentum for tritium atoms ``t``.

The reason for this convention is the existence of the inverse reactions:
``t + d+ -> t+ + d`` outputs diagnostics ``Ftd+_cx`` and ``Fd+t_cx``.

Here are the definitions of some common channels:

+------------------+--------------------------------------+
| Short name       |   Description                        |
+==================+======================================+
| iz               |   Ionisation                         |
+------------------+--------------------------------------+
| ex               |   Excitation (ionisation) radiation  |
+------------------+--------------------------------------+
| rec              |   Recombination                      |
+------------------+--------------------------------------+
| cx               |   Charge exchange                    |
+------------------+--------------------------------------+

.. _sec-flow-diagnostics:

Flow diagnostics convention
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hermes-3 saves a number of different flow rates at a cell boundary. The boundary is 
located at either the "xlow" or "ylow" side (see :ref:`sec-cell-face-notation`).
Positive values correspond to flow into the cell.

There are several types of flow diagnostics at the moment:

+------------------+---------------------------+-------------------------+
| Variable prefix  |   Units                   | Description             |
+==================+===========================+=========================+
| pf               |   :math:`s^{-1}`          | Particle flow           |
+------------------+---------------------------+-------------------------+
| ef               |   :math:`W`               | Energy flow             |
+------------------+---------------------------+-------------------------+
| mf               |   :math:`N`               | Momentum flow           |
+------------------+---------------------------+-------------------------+

There are numerous variants of the above throughout the code, which are
described in their metadata. Note that the diagnostics have not yet been 
implemented for every term.



.. _sec-energy-channels:

Energy channel diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~

Several components can save each term of the energy equations separately, so a simulation
can be power balanced term by term (following Dudson et al. 2026, Nucl. Fusion 66 036015,
section 5 and appendix B). They are named ``E<species>_<term>``, e.g. ``Ed+_conduction``, and
are power densities in :math:`W m^{-3}`. Each one is saved from exactly the field that goes
into the equation, so the channels of an equation add up to its time derivative to rounding.

**Sign.** Positive is a source for whatever the channel belongs to: the species for
``E<species>_*``, and the perpendicular flow energy for ``E<species>_phi_Div_*`` and
``Evort_*``.

**Thermal channels** are :math:`3/2` times the term in the pressure equation, so they are
energy, not pressure. With :math:`\frac{3}{2}` ``SP<species>`` they add up to
:math:`\frac{3}{2}` ``ddt(P<species>)`` (plus ``E<species>_nvh`` when
``numerical_viscous_heating`` is on). In the default form ``E<species>_Div_Vpar_32P`` holds all
of :math:`-\frac{5}{2}\nabla_{||}(P V)` and ``E<species>_VgradP`` holds
:math:`V \partial_{||} P`. With ``p_div_v = true`` they are
:math:`-\frac{3}{2}\nabla_{||}(P V)` and ``E<species>_PdivV``.

**Kinetic channels** are :math:`V` times each term of the parallel momentum equation, so with
:math:`V` ``SNV<species>`` they add up to :math:`V` ``ddt(NV<species>)``. That isn't quite the
rate of change of :math:`\frac{1}{2} m n V^2`, since it leaves out
:math:`-\frac{1}{2} m V^2 \partial n / \partial t`, which matters for electrons (their density
comes from quasineutrality). ``E<species>_Div_VExB_KE_flux`` and ``E<species>_Div_Vpar_KE_flux``
are the conservative flux forms of the two advection terms. They are saved for comparison and
are not part of the sum.

**Volume integrals.** ``E<species>_Div_VExB_32P``, ``E<species>_Div_Vpar_32P``,
``E<species>_Div_Vdia_52P`` (with the default ``diamag_form = 1``) and the two ``_flux``
channels are divergences of a flow. Their volume integral is that flow through the domain
boundary (the targets and the radial boundaries), not a loss inside the domain. The same goes
for each species' pair ``E<species>_Jdia_Grad_phi + E<species>_phi_Div_Jdia``, which is
:math:`-\nabla\cdot(\phi \mathbf{J}_{dia})`. The other kinetic channels and ``Evort_ExB`` are a
velocity or a potential times a divergence, so they only become boundary terms together with
other channels (appendix B).

.. list-table::
   :header-rows: 1
   :widths: 22 26 52

   * - Component
     - Switched on by
     - Channels
   * - ``sheath_boundary_simple``
     - ``diagnose`` in ``[sheath_boundary_simple]``
     - ``E<sp>_sheath_extra`` (see :ref:`sec-sheath_boundary_simple`)
   * - ``evolve_pressure``
     - ``diagnose`` in ``[<sp>]``
     - ``E<sp>_Div_VExB_32P``, ``E<sp>_Div_Vpar_32P``, ``E<sp>_numerical_P`` (all the low
       density, low temperature and low pressure diffusion, hyper diffusion and damping terms)
   * - ``braginskii_conduction``
     - ``diagnose`` in ``[<sp>]``
     - ``E<sp>_conduction``
   * - ``anomalous_diffusion``
     - ``diagnose`` in ``[<sp>]``
     - ``E<sp>_anomalous_D`` (particle diffusion), ``E<sp>_anomalous_chi`` (thermal diffusion)
   * - ``diamagnetic_drift``
     - ``diagnose`` in ``[diamagnetic_drift]``
     - ``E<sp>_Div_Vdia_52P``
   * - ``evolve_momentum``
     - ``diagnose`` in ``[<sp>]``
     - ``E<sp>_Div_VExB_KE``, ``E<sp>_Div_Vpar_KE``, ``E<sp>_V_Gradpar_P``,
       ``E<sp>_j_Gradpar_phi``, ``E<sp>_numerical_V``, and the flux forms
       ``E<sp>_Div_VExB_KE_flux`` and ``E<sp>_Div_Vpar_KE_flux``
   * - ``relax_potential``
     - ``diagnose`` in ``[relax_potential]``
     - ``E<sp>_Jdia_Grad_phi``, ``E<sp>_phi_Div_Jdia``, ``E<sp>_phi_Div_jpar``,
       ``Evort_DivJcol``, ``Evort_ExB``, ``Evort_numerical``, ``Evort_Pi_DivJ``,
       ``Evort_relax`` and ``Vort_mismatch``

**relax_potential.** As in appendix B, the terms of the vorticity equation are weighted by
:math:`-(\phi + \hat{P}_i)`, where :math:`\hat{P}_i` (``Pi_hat`` in the code) is the ion
pressure part of the generalised vorticity. The per species channels take the :math:`-\phi`
part of each current. ``E<sp>_Jdia_Grad_phi`` is the :math:`-\mathbf{J}_{dia}\cdot\nabla\phi`
that goes into the species ``energy_source``. ``Evort_Pi_DivJ`` is the :math:`-\hat{P}_i` part
of the current terms, which nothing in the species equations balances on this path.
``Evort_relax`` is :math:`-(\phi + \hat{P}_i)` times the whole right hand side. With
``evolve_vorticity = false`` it shows how far :math:`\phi` lags behind zero current divergence
(the relaxation residual), and with ``evolve_vorticity = true`` it is the rate of change of the
perpendicular flow energy (up to a divergence). ``Vort_mismatch`` is ``Vort_from_phi - Vort``,
a vorticity rather than a power, and zero when the vorticity isn't evolved. Only the vector
path has channels. On the bracket path (``diamagnetic_bracketform``, the default on FCI meshes)
the component also adds :math:`P \nabla\cdot\mathbf{J}_{dia}` to each species'
``energy_source`` when ``diamagnetic_polarisation`` is set, and that has no channel yet (it
would be ``E<sp>_P_DivJdia``).

**Not covered yet.** A few terms that reach ``energy_source`` or ``ddt`` have no channel, so a
term by term balance leaves them over: the kinetic terms :math:`\frac{1}{2} m R (\Delta v)^2`
of ionisation, recombination and charge exchange, the magnetic flutter terms in
``evolve_pressure`` and ``evolve_momentum``, the electromagnetic :math:`Z A_{||}\partial n/\partial t`
term in ``evolve_momentum``, and ``DivJextra`` in ``relax_potential`` (it only shows up in
``Evort_relax``). ``neutral_mixed`` has no channels.
