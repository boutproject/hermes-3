# Continuum ion energy MMS

BOUT++ expression sources evolve two ions, zero-current electrons, vorticity
and total ion energy on a physical-X, periodic-Y/Z curvature slab. Baseline,
skew metric (`g13=1/4`, `g33=17/16`) and viscous (`nu=0.02`) variants run at
sizes 8 and 16 with twice as many Y cells, through normalised time 0.1. Flat parallel-velocity
extrema avoid loss of order in limited transport. Analytic curvature-driven
polarisation flow exercises every evolved variable; the continuum phi is zero
and its recovery is checked. Every variable and recovered phi, pressure
potential and ion temperature must have L2 order at least 1.7 and satisfy the
finest-grid error limits. A repeated baseline case checks diagnostics on/off.
The `long` profile adds 32 and requires order 1.8 on the finest pair. CTest runs
it on two MPI ranks and reserves two processor slots. The skew case uses a
Dirichlet auxiliary AC boundary; its analytic value and wall flux both vanish.

Known failure: `ion-energy-mms-long` fails the skew vorticity 16-to-32 check:
measured L2 order 1.725, required 1.8. All other long-profile checks pass,
including baseline and viscosity convergence. The default profile passes.
The failing check remains enabled with its original threshold; its investigation
is deferred.

```sh
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 --profile long --nproc 2
python tests/integrated/ion-energy-mms/mms.py --all --check --check-files
```

The runner uses committed inputs and needs neither SymPy nor plotting packages.
SymPy is required only to regenerate inputs and verify continuum identities.
The default profile took 216 seconds in the fresh-build validation with two
concurrent CTest jobs and one MPI rank per case.
Recovered pressure uses mean-Neumann, zero-fluctuation midpoint wall guards.
Parallel current, varying magnetic field and Reynolds advection are not covered
by these profiles. The closed test checks nontrivial perpendicular transport.
Results include orders, errors, recovery, hashes and execution times in
`result.json`.
