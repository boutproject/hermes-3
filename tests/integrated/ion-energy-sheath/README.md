# Production sheath and parallel-current regression

Run `runtest --nproc 2` with the configured energy-analysis Python environment.
The small all-open identity slab has physical targets at both Y ends; two
ranks split Y with `NXPE=1`. It evolves electron momentum independently of ion
momentum and advances to normalised time 0.02. This is a boundary/current
regression, not a turbulence or closed-energy test.

The test checks target currents against independent Bohm/electron sheath
formulas at the initial state, reconstructs current divergence and completed
parallel electrical work from saved completed-current and potential guards, and checks
auxiliary continuity at every physical cell including target-adjacent cells.
The two-rank profile compares all physical histories with a serial run to
`1e-10` absolute error.
Limits are declared in `criteria.json`; no existing MMS accuracy gate applies.

The companion C++ production-sheath test perturbs only electron momentum
guards and verifies the analytic current-divergence and auxiliary response,
including extra/collisional currents. Physical X walls allow the net parallel
current source to close; periodic-X unsolvable means are tested separately.

The maintainer accepts the previously observed component-reordering differences
of up to `8.881784197001252e-16` as roundoff. The redundant reordered run and
bitwise sorting assertions have been removed. Target, continuity, electrical
and Y-MPI checks retain their thresholds. Each invocation saves its measured
checks and serial/parallel execution details in `result.json`.
