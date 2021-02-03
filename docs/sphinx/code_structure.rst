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

    * `collision_frequency`   Normalised collision frequency
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
  
Docs
----

The basic building block of all Hermes-3 models is the
`Component`. This defines an interface to a class which takes a state
(a tree of dictionaries/maps), and transforms (modifies) it.  After
all components have modified the state in turn, all components may
then implement a `finally` method to take the final state but not
modify it. This allows two components to depend on each other, but
makes debugging and testing easier by limiting the places where the
state can be modified.

.. doxygenstruct:: Component
   :members:

Components are usually defined in separate files; sometimes multiple
components in one file if they are small and related to each other (e.g.
atomic rates for the same species). To be able to create components,
they need to be registered in the factory. This is done in the header
file using a code like::

  #include "component.hxx"

  struct MyComponent : public Component {
    ...
  };
  
  namespace {
  RegisterComponent<MyComponent> registercomponentmine("mycomponent");
  }

where `MyComponent` is the component class, and "mycomponent" is the
name that can be used in the BOUT.inp settings file to create a
component of this type. Note that the name can be any string except it
can't contain commas or brackets (), and shouldn't start or end with
whitespace.
