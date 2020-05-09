.. _sec-code_structure:

Code structure
==============

Component initialisation
------------------------

Inputs to the constructors are:

* `name`
* `alloptions`
* `solver`

The `name` is a string labelling the instance. The `alloptions` tree contains at least:

* `alloptions[name]` options for this instance
* `alloptions['units']` 

State
-----

The simulation state is passed between components, and is
a tree of objects (Options). At the start of each iteration
(rhs call) it contains:

* `time`   BoutReal, the current simulation time
* `units`
  
  * `seconds`   Multiply by this to get units of seconds
  * `eV`          Temperature normalisation
  * `Tesla`       Magnetic field normalisation
  * `meters`      Length normalisation
  * `inv_meters_cubed`     Density normalisation

* `species`  Plasma species

  * `e`
  * `species1`  Example "h", "he2+"

    * `AA`  Atomic mass, proton = 1
    * `charge`  Charge, in units of proton charge (i.e. electron=-1)
    
    * `density`
    * `momentum`
    * `pressure`
    * `velocity` Parallel velocity
    * `temperature`

    * `collision_rate`   Normalised collision frequency
    * `density_source`
    * `momentum_source`
    * `energy_source`

* `fields`

  * `vorticity`
  * `phi`       Electrostatic potential
  * `DivJdia`   Divergence of diamagnetic current
  * `DivJextra` Divergence of current, including parallel current.
    Not including diamagnetic or polarisation currents


Notes:
* When checking if a subsection exists, use `option.isSection`, since `option.isSet`
  is false if it is a section and not a value.
  
