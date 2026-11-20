# Closed two-ion CI regression

An 8 x 4 x 8 periodic slab evolves two ion species with varying composition,
parallel momenta, curvature and viscosity. Fields are Y-invariant: the
zero-current electron closure does not provide general parallel force balance.
Parallel kinetic energy changes through perpendicular transport.

The two diagnostics-on/off runs check plasma energy, species particles and
charge, energy recovery, reconstructed diagnostics and state agreement at the
unchanged `1e-10` limit. Curvature and perpendicular kinetic exchanges must each
exceed `1e-4` of initial energy; parallel kinetic exchange must exceed `5e-5`.
The run ends at t=20; reducing resolution preserves all acceptance criteria.

```sh
ctest --test-dir build -R '^ion-energy-closed$' --output-on-failure
python tests/integrated/ion-energy-closed/runtest --exe build/hermes-3
# Optional resolution check
python tests/integrated/ion-energy-closed/runtest --exe build/hermes-3 --resolution 16
```

The default uses two MPI ranks. An 8 x 4 x 8 benchmark completed both runs and
analysis in 13.2 s, with relative energy drift 5.5e-16. CTest reserves two cores
and has a 60 s ceiling. `result.json` records measured checks, hashes and timing.
There are no sources, sheaths or parallel current. Source conversion, sheath
equations and boundary identities retain production unit-test coverage.
