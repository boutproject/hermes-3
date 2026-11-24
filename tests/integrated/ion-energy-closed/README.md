# Closed ion-energy regression

A periodic, Y-invariant two-ion slab checks conservation of `Eion + 1.5 Pe`
and each ion inventory to relative 1e-10, nonzero thermal/kinetic exchange,
thermodynamic recovery from saved production diagnostics, and identical states
with diagnostics enabled or disabled. The `zero_current` electron closure has
no independently evolved kinetic reservoir in this configuration.

```sh
python tests/integrated/ion-energy-closed/runtest --exe build/hermes-3 \
  --output-root build/closed-check
```

Use a fresh output directory. The default is two ranks and 8x4x8 physical
cells; `--resolution` changes X and Z. Requires boutdata and xhermes on its
`total-energy` branch. CTest runs this regression with a 60 second timeout.
