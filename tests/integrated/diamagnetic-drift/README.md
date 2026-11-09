# Diamagnetic drift on a single-null tokamak

This axisymmetric test uses the repository's COMPASS 36x48 grid, metric and
curvature. A small BOUT++ model evolves density, pressure and parallel momentum
using only the public `diamagnetic_drift` component. Other transport, sources
and sheath losses are absent. Temperature is P/N, with deliberately different
target guard values to exercise sheath-face closure.
The driver communicates the grid's magnetic-field guards before constructing
the component, so normalised curvature agrees across processor interfaces.

The divergence form with `closed_average` is compared with `closed` and an
independent finite-volume reference that explicitly makes every physical face
velocity zero, then calls the operator with radial flux enabled. The comparison
checks every cell for all three sources. The expected core average is computed
separately; private-flux cells must retain their unaveraged source. Radial guard
fluxes must vanish pointwise. Positive controls require nonzero flux through
the core, private-flux wall, outer wall and both targets if those faces are open.
These checks prevent opposing boundary leaks from cancelling in a global sum.

Every output checks volume-weighted density, pressure and momentum conservation
and zero source integrals to a relative tolerance of 1e-12. The core ring stays
uniform. Fields and source diagnostics agree across processor layouts to a
relative maximum-norm tolerance of 1e-10. Evolution must be nonzero.

COMPASS has eight cells in each divertor leg, so BOUT++ requires at least six
Y processors to put its branch cuts on processor boundaries. The test compares
a six-rank reference (`NXPE=1, NYPE=6`) with twelve-rank layouts `(1,12)` and
`(2,6)`, rather than attempting invalid one-, two- or four-rank runs. Each run
has a 60-second timeout that also terminates its MPI process group.

Build `diamagnetic_drift_smoke`, then run `ctest -R '^diamagnetic-drift$'` from
the build directory with the project's Python environment on `PATH`.
CTest passes the executable, grid and output paths explicitly, including when
`HERMES_COPY_TESTS_TO_BUILD` is disabled. The script also accepts `--executable`,
`--grid`, `--run-root` and `--timeout`. It uses `MPIRUN` if set, otherwise
`mpiexec -np`.

The real grid's curvature is used with the conservative divergence form.
Gradient-form conservation is covered by unit tests with a curvature whose
discrete divergence vanishes for C2 derivatives.
