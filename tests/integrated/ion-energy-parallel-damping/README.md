# Parallel potential damping probe

Run `runtest --exe /path/to/hermes-3 --output-root /path/to/evidence` with the
configured Python environment. This small unforced periodic slab evolves the
real ion-energy, conservative drift and electron-momentum components for 0.02.
An alternating-Y potential with a nonzero Z mode has zero collocated C2
parallel derivative and a nonzero numerical face damping rate.

The runner compares damping off/on at an identical initial state, checks the
completed physical rate and auxiliary continuity, records evolved ion energy
and its recovered thermal/perpendicular partition, and compares serial and
two-rank Y decomposition with the existing `1e-10` tolerance. Numerical damping
uses collective sound speed and contributes no explicit energy source. The
production guards and equations are used throughout.
