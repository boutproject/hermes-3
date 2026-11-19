# Ion-energy radial-wall regressions

`ion-energy-frozen/BOUT.inp` defines an exact discrete state using input
expressions. Its conducting potential is
φ = 0.01 sin(πx) cos(z) + 0.01 sin(2πx), and its Neumann pressure is
p = ψ = 1 + 0.05 cos(πx) sin(z). Vorticity is the finite-volume Laplacian
of Φ = φ + ψ, and Eion = 3p/2 + K⊥ uses the same face differences.
The test doubles 16² to 32² and 64² on one rank. For both walls, adjacent
differences of the native ddt(Eion) and DivJReynolds RMS values must shrink
by a factor between 3.5 and 4.5. The relative potential and absolute vorticity
recovery tolerances are 1e-12, with potential denominator scale 1e-3 eV.

`ion-energy-multispecies/BOUT.inp` is a 16×8×16 periodic-Y slab with two
Y guards. Deuterium and helium have varying wall fractions and nonzero,
varying parallel velocities. The input expressions define consistent discrete
vorticity, perpendicular energy and parallel energy. Through t = 0.1 the
test checks each species' particle and momentum face-source identity, the
energy/current identities, local thermal and kinetic weights, zero unpaid
wall work and positive ion temperature.

The slow tests generate the fixed-wall scan input from
`examples/tokamak-2D-driftplane/2D-energy/scan-input.inp.template`, using
its matched-startup generator and the original solver tolerances. Both use
custom_D = 0.1 and viscosity = 1. The convergence test uses two MPI ranks,
matching the example comparison; the wall-envelope test uses one rank.
The pressure perturbation has
fixed physical wall buffers and a C² taper, so refinement samples the same
physical initial condition. Generated restarts are runtime artifacts.

The convergence test runs 16²–64² through t = 5000 at output spacing 100.
The second adjacent K⊥ difference must be less than half the first; the
outer-10% kinetic-energy share must vary by less than 0.005 and grow by no
more than 0.001 from the coarsest to finest grid. Exact partial-cell overlaps
define the strip at every resolution.

The wall-amplitude test runs 128² through t = 1200 at spacing 5. Its
`ion-energy-wall-envelope.json` records the generated input hash, resolved
configuration, measured maxima and bounds. Native version-1 maxima are
2.11665e-8 in modes 16–31 and 3.13835e-8 in modes 32–63. Bounds of 3.5e-8
and 5e-8 give about 60% margin. These are absolute first-cell wall RMS values,
independent of interior amplitude. The input deliberately matches the
example physical setup, including classical diffusion, rather than the
small undiffused gauge input.
The envelope retains the measured input hash separately: renaming the runner
changes only the current input's header comment, with identical physical options.

Run the native checks after building Hermes:

```sh
PYTHONPATH=tools:external/BOUT-dev/tools/pylib python -m pytest -q tests/python/test_local_wall_transport.py tests/python/test_ion_energy_native.py
PYTHONPATH=tools:external/BOUT-dev/tools/pylib python -m pytest -q --run-slow tests/python/test_ion_energy_convergence.py tests/python/test_ion_energy_native.py -m slow
```

Set `HERMES_EXECUTABLE` to select another executable and `--basetemp` to retain
native runs in a chosen directory. The slow tests write `convergence.json`
and `wall-bands.json` beside their native outputs. Regenerate the envelope
from a fresh complete wall trajectory; retain its configuration and measured
maxima when changing bounds.
