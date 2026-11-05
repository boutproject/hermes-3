# Quadratic wall relaxation

Run `./runtest --exe /path/to/hermes-3` with the project Python dependencies.
Use `--output-root /path/to/results` to choose where case directories are saved.

This test reuses the orthogonal slab input from the mixed-polarisation-boundary
test, disables curvature, and starts uniform positive or negative vorticity with
zero wall midpoint values. It runs four output intervals at 16, 32 and 64
interior X cells over the same physical width.

Both walls must have the analytic quadratic drift rate
`-kappa_hat * sqrt(g11) * omega * Lx / (2 * c * g11)`. The saved wall RHS, the
change in saved phi midpoints between outputs, and the auxiliary chi must agree
with that rate to relative tolerance `1e-6`. The small vorticity amplitude keeps
pressure-potential feedback below this budget. The former one-sided stencil
fails this test with errors of 4.17%, 2.08% and 1.04%.
