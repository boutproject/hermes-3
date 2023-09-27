.. _sec-tokamak_turbulence

Tokamak turbulence
==================

Turbulence of turbulence in 3D tokamak geometries. Hermes-3 solves an
electrostatic 6-field model for vorticity, electron density, electron
and ion parallel velocity, electron and ion pressure.

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
density by quasineutrality.

   
Preparing an input mesh
-----------------------

Generate a mesh file using Hypnotoad

Adjusting curvature



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
