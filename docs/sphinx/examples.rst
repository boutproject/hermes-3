.. _sec-examples:

Examples
========


1D flux-tube
------------

These simulations follow the dynamics of one or more species along the
magnetic field. By putting a source at one end of the domain, and a
sheath at the other, this can be a useful model of plasma dynamics in
the Scrape-Off Layer (SOL) of a tokamak or other magnetised plasma.

1D periodic domain, Te and Ti
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A fluid is evolved in 1D, imposing quasineutrality and zero net current.
Both electron and ion pressures are evolved, but there is no exchange
of energy between them, or heat conduction.

.. figure:: figs/1d_te_ti.*
   :name: 1d_te_ti
   :alt:
   :width: 60%
   
   Evolution of pressure, starting from a top hat

The model components are ions (i), electrons (e), and a constraint
that the net current is zero. This constraint applies the electron
pressure to the ion momentum equation, and sets the electron velocity
to be equal to the ion velocity.

.. code-block:: ini

   [hermes]
   components = i, e, zero_current


The ion density, pressure and momentum equations are evolved:

.. code-block:: ini

   [i]  # Ions
   type = evolve_density, evolve_pressure, evolve_momentum


The electron density is set to the ion density by quasineutrality,
and only the electron pressure is evolved.

.. code-block:: ini

   [e] # Electrons
   type = quasineutral, evolve_pressure


2D drift-plane
--------------

Simulations where the dynamics along the magnetic field is not
included, or only included in a parameterised way as sources or
sinks. These are useful for the study of the basic physics of plasma
"blobs" / filaments, and tokamak edge turbulence.

Blob2d
~~~~~~

A seeded plasma filament in 2D. This version is isothermal and cold ion,
so only the electron density and vorticity are evolved. A sheath-connected
closure is used for the parallel current.

.. figure:: figs/blob2d.png
   :name: blob2d
   :alt:
   :scale: 50
   
   Electron density Ne at three times, showing propagation to the right

The model components are

.. code-block:: ini

   [hermes]
   components = e, vorticity, sheath_closure

The electron component consists of two types:

.. code-block:: ini

   [e]  # Electrons
   type = evolve_ne, isothermal


The `evolve_ne` component type evolves the electron density `Ne`. This component
has several options, which are set in the same section e.g.

.. code-block:: ini

   poloidal_flows = false  # Y flows due to ExB

The `isothermal` component type sets the temperature to be a constant, and using
the density then sets the pressure. The constant temperature is also
set in this `[e]` section:

.. code-block:: ini

   temperature = 5  # Temperature in eV

The `vorticity` component uses the pressure to calculate the diamagnetic current,
so must come after the `e` component. This component then calculates the potential.
Options to control the vorticity component are set in the `[vorticity]` section.

The `sheath_closure` component uses the potential, so must come after `vorticity`.
Options are also set as

.. code-block:: ini

   [sheath_closure]
   connection_length = 10 # meters

Blob2D-Te-Ti
------------

A seeded plasma filament in 2D. This version evolves both electron and
ion temperatures. A sheath-connected closure is used for the parallel
current.

.. figure:: figs/blob2d-te-ti.png
   :name: blob2d-te-ti
   :alt:
   :scale: 50
   
   Electron density Ne at three times, showing propagation to the right and downwards

The model components are

.. code-block:: ini

   [hermes]
   components = e, h+, vorticity, sheath_closure


The electron component evolves density (saved as `Ne`) and pressure
(`Pe`), and from these the temperature is calculated.

.. code-block:: ini

   [e]
   type = evolve_density, evolve_pressure


The ion component sets the ion density from the electron density, by
using the quasineutrality of the plasma; the ion pressure (`Ph+`) is evolved.

.. code-block:: ini
   
   [h+]
   type = quasineutral, evolve_pressure


2D axisymmetric tokamak
-----------------------

These are transport simulations, where the cross-field transport is given
by diffusion, and fluid-like equations are used for the parallel dynamics
(as in the 1D flux tube cases).

The `recycling-dthene` example includes cross-field diffusion,
parallel flow and heat conduction, collisions between species, sheath
boundary conditions and recycling. It simulates the density, parallel
flow and pressure of the electrons; ion species D+, T+, He+, Ne+; and
neutral species D, T, He, Ne.

.. figure:: figs/pe_nvt_nne_2d.png
   :name: recycling-dthene
   :alt:
   :width: 100%

   Electron pressure, parallel tritium flux, and neon atom density. Simulation
   evolves D, T, He, Ne and electron species, including ions and neutral atoms.

The model components are a list of species, and then collective components
which couple multiple species.

.. code-block:: ini

   [hermes]
   components = (d+, d, t+, t, he+, he, ne+, ne, e,
                 collisions, sheath_boundary, recycling, reactions)

Note that long lists like this can be split across multiple lines by
using parentheses. 
                 
Each ion species has a set of components, to evolve the density,
momentum and pressure. Anomalous diffusion adds diffusion of
particles, momentum and energy. For example deuterium ions contain:

.. code-block:: ini
   
   [d+]
   type = evolve_density, evolve_momentum, evolve_pressure, anomalous_diffusion
   AA = 2
   charge = 1

Atomic reactions are specified as a list:

.. code-block:: ini
   
   [reactions]
   type = (
        d + e -> d+ + 2e,   # Deuterium ionisation
        t + e -> t+ + 2e,   # Tritium ionisation
        he + e -> he+ + 2e, # Helium ionisation
        he+ + e -> he,      # Helium+ recombination
        ne + e -> ne+ + 2e, # Neon ionisation
        ne+ + e -> ne,      # Neon+ recombination
       )
