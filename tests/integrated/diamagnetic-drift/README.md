# Diamagnetic particle conservation on a tokamak grid

Run Hermes-3 with `diamagnetic_drift` and a single species evolving density
at prescribed temperature. There is no parallel velocity, sheath loss or
external source. The bundled `examples/tokamak-3D/tokamak.nc` supplies the
geometry and curvature; no transport operator is copied into the test.

For MZ=1 and MZ=4, `closed` and `closed_average` must conserve the
volume-weighted particle content to relative 1e-12. `open` must change it,
and every case must have nonzero density evolution. Physical X and Y guards
are excluded from the integrals. Six Y ranks exercise the core averaging
across MPI processes and place the grid's divertor branch cuts correctly.

```sh
python tests/integrated/diamagnetic-drift/runtest --exe build/hermes-3 \
  --grid examples/tokamak-3D/tokamak.nc --output-root build/diamagnetic-check
ctest --test-dir build -R '^diamagnetic-drift$' --output-on-failure
```

Requires boutdata. Launches use `launch_safe`; CTest supplies executable and
grid paths and a timeout. `--output-root` selects a fresh directory; otherwise
inputs, outputs and logs are retained in a temporary directory.
