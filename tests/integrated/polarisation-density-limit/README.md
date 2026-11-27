# Polarisation density transport regression

The existing closed two-ion case checks that zero and negligible density
thresholds produce the same state to absolute 1e-10. Supply saved base-branch
output with `--baseline-directory` to compare against that branch as well.

The wall case uses an input-expression density band at the inner radial wall.
It checks active limiter diagnostics, particle and energy wall flux closure,
the saved total-ion-energy RHS balance, positive density over a fixed run,
and the final inner-wall density with the limiter off and on. It also checks
optional output fields, central and upwind weights, and independence of radial
MPI decomposition and diagnostic output settings.

```sh
PYTHONPATH=external/xhermes python tests/integrated/polarisation-density-limit/runtest \
  --exe build/hermes-3 --output-root build/polarisation-runs
```

Use a fresh output directory. Requires boutdata and xhermes with its
total-energy inventory API. The default is two ranks; `--nproc` controls the
primary decomposition. Native outputs and launch logs are retained in the
output directory. Unit tests separately check face-rate values with changing
composition, both current signs, and core/PF classification.
