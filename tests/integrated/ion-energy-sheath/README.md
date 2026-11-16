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
Reversed component construction must agree bitwise. The two-rank profile also
compares all physical histories with a serial run to `1e-10` absolute error.
Limits are declared in `criteria.json`; no existing MMS accuracy gate applies.

The companion C++ production-sheath test perturbs only electron momentum
guards and verifies the analytic current-divergence and auxiliary response,
including extra/collisional currents. Physical X walls allow the net parallel
current source to close; periodic-X unsolvable means are tested separately.

Current validation: the maintainer accepts the observed reversed-construction
history difference `8.881784197001252e-16` as probably harmless ordering roundoff.
The original zero thresholds and 14 raw bitwise failures remain retained; the
acceptance covers those measured histories only. Identical-state
C++ RHS comparisons pass bitwise; all target, continuity, electrical and Y-MPI
checks pass. See `docs/polarisation-drift-3d-analysis-report.md` for retained raw
results. Thresholds remain unchanged.
