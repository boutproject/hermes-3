.. _sec-code_structure:

Code structure
==============



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


