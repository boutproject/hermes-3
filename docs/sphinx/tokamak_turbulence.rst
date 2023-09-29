.. _sec-tokamak_turbulence

Tokamak turbulence
==================

Turbulence of turbulence in 3D tokamak geometries. Hermes-3 is
configured to solve an electrostatic 6-field model for vorticity,
electron density, electron and ion parallel velocity, electron and ion
pressure.

The input file is in the Hermes-3 repository under
``examples/tokamak/turbulence``.

Model equations
---------------

The lines that define the components to include in the model are:

.. code-block:: ini

   [hermes]
   components = (e, d+, sound_speed, vorticity,
                 sheath_boundary, collisions,
                 diamagnetic_drift, classical_diffusion,
                 polarisation_drift
                )

   [e]
   type = evolve_density, evolve_pressure, evolve_momentum

   [d+]
   type = quasineutral, evolve_pressure, evolve_momentum

We define two species: electrons ``e`` and deuterium ions ``d+``.
Electron density is evolved, and ion density is set to electron
density by quasineutrality.  The electron fluid equations for density
:math:`n_e`, parallel momentum :math:`m_en_ev_{||e}`, and pressure
:math:`p_e = en_eT_e` are:

.. math::

   \begin{aligned}
   \frac{\partial n_e}{\partial t} =& -\nabla\cdot\left[n_e \left(\mathbf{v}_{E\times B} + \mathbf{b}v_{||e} + \mathbf{v}_{de}\right)\right] + S_n \\
   \frac{\partial}{\partial t}\left(m_en_ev_{||e}\right) =& -\nabla\cdot\left[m_en_ev_{||e} \left(\mathbf{v}_{E\times B} + \mathbf{b}v_{||e} + \mathbf{v}_{de}\right)\right] - \mathbf{b}\cdot\nabla p_e \nonumber \\
   &- en_eE_{||} + F_{ei} \\
   \frac{\partial}{\partial t}\left(\frac{3}{2}p_e\right) =& -\nabla\cdot\left[\frac{3}{2}p_e \left(\mathbf{v}_{E\times B} + \mathbf{b}v_{||e}\right) + \frac{5}{2}p_e\mathbf{v}_{de}\right] - p_e\nabla\cdot\left(\mathbf{v}_{E\times B} + \mathbf{b}v_{||e}\right) \nonumber \\
   & + \nabla\cdot\left(\kappa_{||e}\mathbf{b}\mathbf{b}\cdot\nabla T_e\right) + S_{Ee} + W_{ei}
   \end{aligned}

Here the electrostatic approximation is made, so :math:`E_{||} = -\mathbf{b}\cdot\nabla\phi`.

The ion fluid equations assume quasineutrality so :math:`n_i = n_e`,
and evolve the ion parallel momentum :math:`m_in_iv_{||i}` and
pressure :math:`p_i`:

.. math::

   \begin{aligned}
   \frac{\partial}{\partial t}\left(m_in_iv_{||i}\right) =& -\nabla\cdot\left[m_in_iv_{||i} \left(\mathbf{v}_{E\times B} + \mathbf{b}v_{||i} + \mathbf{v}_{di}\right)\right] - \mathbf{b}\cdot\nabla p_i \nonumber \\
   &+ Z_ien_iE_{||} - F_{ei} \\
   \frac{\partial}{\partial t}\left(\frac{3}{2}p_i\right) =& -\nabla\cdot\left[\frac{3}{2}p_i \left(\mathbf{v}_{E\times B} + \mathbf{b}v_{||i}\right) + \frac{5}{2}p_i\mathbf{v}_{di}\right] - p_i\nabla\cdot\left(\mathbf{v}_{E\times B} + \mathbf{b}v_{||i}\right) \nonumber \\
   & + \nabla\cdot\left(\kappa_{||i}\mathbf{b}\mathbf{b}\cdot\nabla T_i\right) + S_{Ei} + S_n\frac{1}{2}m_in_iv_{||i}^2 - W_{ei} \nonumber \\
   & + \frac{p_i}{en_0}\nabla\cdot\left(\mathbf{J}_{||} + \mathbf{J}_d\right)
   \end{aligned}

The vorticity :math:`\omega` is

.. math::

   \omega = \nabla\cdot\left[\frac{m_in_0}{B^2}\nabla_\perp\left(\phi + \frac{p_i}{n_0}\right)\right]

whose evolution is given by the current continuity equation:

.. math::

   \begin{aligned}
   \frac{\partial \omega}{\partial t} =& -\nabla\cdot\left[\frac{m_i}{2B^2}\nabla_\perp\left(\mathbf{v}_E \cdot\nabla p_i\right) + \frac{\omega}{2}\mathbf{v}_E + \frac{m_in_0}{2B^2}\nabla_\perp^2\phi\left(\mathbf{v}_E + \frac{\mathbf{b}}{n_0B}\times\nabla p_i\right)\right] \nonumber \\
   &+ \nabla\cdot\left(\mathbf{J}_{||} + \mathbf{J}_d + \mathbf{J}_{ci}\right)
   \end{aligned}

where the Boussinesq approximation is made, replacing the density in
the polarisation current with a constant :math:`\overline{n}`.  The
divergence of the diamagnetic current is written as

.. math::

   \nabla\cdot\mathbf{J}_d = \nabla\cdot\left[\left(p_e + p_i\right)\nabla\times\frac{\mathbf{b}}{B}\right]
   
Preparing an input mesh
-----------------------

Generate a mesh file using Hypnotoad

Adjusting curvature


Compiling Hermes-3
------------------



Starting a simulation
---------------------

The input file sets the number of output steps to take and the time
between outputs in units of reference ion cyclotron time:

.. code-block:: ini

   nout = 10      # Number of output steps
   timestep = 10  # Output timestep, normalised ion cyclotron times [1/Omega_ci]

With the normalisations in the input that reference time is about 1e-8
seconds, so taking 10 steps of 10 reference cyclotron times each
advances the simulation by around 1 microsecond in total.

The first few steps are likely to be slow, but the simulation should
speed up considerably by the end of these 10 steps. This is largely
due to rapid transients as the electric field is set up by the sheath
and parallel electron flows.

Adjusting input sources
-----------------------

**Note**: When starting a new simulation, it is important to calibrate
the input sources, to ensure that the particle and power fluxes are
what you intend.

The inputs are electron density source, electron and ion heating power.
The particle source is set in the electron density section ``[Ne]``:

.. code-block:: ini

   [Ne]
   flux = 3e21 # /s
   shape_factor = 1.0061015504152746

   source = flux * shape_factor * exp(-((x - 0.05)/0.05)^2)
   source_only_in_core = true

The inputs read by Hermes-3/BOUT++ are ``source`` and
``source_only_in_core``. The ``flux`` and ``shape_factor`` values are
just convenient ways to calculate the source (New variables can be
defined and used, and their order in the input file doesn't matter).

Continuing the simulations
--------------------------

A turbulence simulation typically takes many days of running, to reach
(quasi-)steady state then gather statistics for analysis.  To continue
a simulation, the simulation state is loaded from restart files
(BOUT.restart.*) and the simulation continues running. The "nout" and
"timestep" set the number of *new* steps to take. To do this, copy
the BOUT.inp (options) file and BOUT.restart.* files into a new directory.
For example, if the first simulation was in a directory "01":

.. code-block:: bash

   $ mkdir 02
   $ cp 01/BOUT.inp 02/
   $ cp 01/BOUT.restart.* 02/

We now have a new input file (02/BOUT.inp) that we can edit to update
settings. I recommend increasing the output ``timestep`` from 10 to 100,
and the number of outputs ``nout`` from 10 to 100. You can also adjust
particle and power sources, or make other changes to the settings. Once
ready, restart the simulation:

.. code-block:: bash

   $ mpirun -np 64 ./hermes-3 -d 02 restart

Note the ``restart`` argument.



