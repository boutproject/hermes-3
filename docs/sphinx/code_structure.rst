.. _sec-code_structure:

Code structure
==============

A hermes-3 model, like all `BOUT++ models
<https://bout-dev.readthedocs.io/en/latest/user_docs/physics_models.htmlject.github.io/>`_,
is an implementation of a set of Ordinary Differential Equations
(ODEs). The time integration solver drives the simulation, calling the
`Hermes::rhs` function to calculate the time-derivatives of all the
evolving variables.

The calculation of the time derivatives is coordinated by passing
a state object between components. The state is a nested tree, and
can have values inserted and retrieved by the components. The components
are created and then run by a scheduler, based on settings in the
input (BOUT.inp) file.

In terms of design patterns, the method used here is essentially a combination
of the `Encapsulate Context <https://accu.org/journals/overload/12/63/kelly_246/>`_
and `Command <https://en.wikipedia.org/wiki/Command_pattern>`_ patterns.

Simulation state
----------------

The simulation state is passed between components, and is a tree of
objects (Options objects). At the start of each iteration (rhs call) a
new state is created and contains:

* `time`   BoutReal, the current simulation time
* `units`
  
  * `seconds`   Multiply by this to get units of seconds
  * `eV`          Temperature normalisation
  * `Tesla`       Magnetic field normalisation
  * `meters`      Length normalisation
  * `inv_meters_cubed`     Density normalisation

so the temperature normalisation can be extracted using::

  BoutReal Tnorm = state["units"]["eV"];
    
As the components of a model are run, they set, modify and use values
stored in this state. To ensure that components use consistent names
for their input and output variables, a set of conventions are used
for new variables which are added to the state:

* `species`  Plasma species

  * `e`    Electron species
  * `species1`  Example "h", "he+2"

    * `AA`  Atomic mass, proton = 1
    * `charge`  Charge, in units of proton charge (i.e. electron=-1)
    
    * `density`
    * `momentum` Parallel momentum
    * `pressure`
    * `velocity` Parallel velocity
    * `temperature`

    * `collision_frequency`   Normalised collision frequency
    * `density_source`  Normalised particle source 
    * `momentum_source` Normalised momentum source
    * `energy_source`  Normalised energy source

    * `particle_flow_xlow` Normalised particle flow through lower X cell face
    * `particle_flow_ylow` Normalised particle flow through lower Y cell face
    * `momentum_flow_xlow` Normalised momentum flow through lower X cell face
    * `momentum_flow_ylow` Normalised momentum flow through lower Y cell face
    * `energy_flow_xlow`   Normalised energy flow through lower X cell face
    * `energy_flow_ylow`   Normalised energy flow through lower Y cell face

* `fields`

  * `vorticity`
  * `phi`           Electrostatic potential
  * `Apar`          Electromagnetic potential b dot A in induction terms
  * `Apar_flutter`  The electromagnetic potential (b dot A) in flutter terms
  * `DivJdia`       Divergence of diamagnetic current
  * `DivJcol`       Divergence of collisional current
  * `DivJextra`     Divergence of current, including 2D parallel current
                    closures.  Not including diamagnetic, parallel current due to
                    flows, or polarisation currents

For example to get the electron density::

  Field3D ne = state["species"]["e"]["density"];

This way of extracting values from the state will print the value to
the log file, and is intended mainly for initialisation. In
`Component::transform` and `Component::finally` functions which run
frequently, faster access methods are used which don't print to the
log. To get a value::

  Field3D ne = get<Field3D>(state["species"]["e"]["density"]);

If the value isn't set, or can't be converted to the given type,
then a `BoutException` will be thrown.

To set a value in the state, there is the `set` function::

  set(state["species"]["h"]["density"], ne);

A common need is to add or subtract values from fields, such as density sources::

  add(state["species"]["h"]["density_source"], recombination_rate);
  subtract(state["species"]["h+"]["density_source"], recombination_rate);
  
Notes:

- When checking if a subsection exists, use `option.isSection`, since `option.isSet`
  is false if it is a section and not a value.
- The species name convention is that the charge state is last, after the `+` or `-`
  sign: `n2+` is a singly charged nitrogen molecule, while `n+2` is a +2 charged
  nitrogen atom.
  
Components
----------

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
    MyComponent(const std::string &name, Options &options, Solver *solver);
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

Inputs to the component constructors are:

* `name`
* `alloptions`
* `solver`

The `name` is a string labelling the instance. The `alloptions` tree contains at least:

* `alloptions[name]` options for this instance
* `alloptions['units']`
  

Component scheduler
-------------------

The simulation model is created in `Hermes::init` by a call to the `ComponentScheduler`::

  scheduler = ComponentScheduler::create(options, Options::root(), solver);

and then in `Hermes::rhs` the components are run by a call::

  scheduler->transform(state);

The call to `ComponentScheduler::create` treats the "components"
option as a comma-separated list of names. The order of the components
is the order that they are run in. For each name in the list, the
scheduler looks up the options under the section of that name. 

.. code-block:: ini

   [hermes]
   components = component1, component2

   [component1]

   # options to control component1

   [component2]

   # options to control component2

This would create two `Component` objects, of type `component1` and
`component2`. Each time `Hermes::rhs` is run, the `transform`
functions of `component1` amd then `component2` will be called,
followed by their `finally` functions.

It is often useful to group components together, for example to
define the governing equations for different species. A `type` setting
in the option section overrides the name of the section, and can be another list
of components

.. code-block:: ini

   [hermes]
   components = group1, component3

   [group1]
   type = component1, component2
   
   # options to control component1 and component2

   [component3]

   # options to control component3

This will create three components, which will be run in the order
`component1`, `component2`, `component3`: First all the components
in `group1`, and then `component3`. 

.. doxygenclass:: ComponentScheduler
   :members:
