# Drift-plane energy publication scan

The scan compares matched legacy pressure and total ion energy formulations.
Both include `classical_diffusion`, with a constant normalised diffusivity
`custom_D=1.19` (17 m²/s for this setup), configurable with `--custom-D`.
It generates refinement tables and figures and supports preparing and importing
remote runs. See [publication-study.md](publication-study.md) for the full study.

`ion-thermal-gain-<boundary>.pdf` plots the signed full-simulation gain
`U_i(T)-U_i(0)` against mesh spacing on a logarithmic axis. Available pressure
and ion-energy runs share the plot; with at least three suitable ion-energy
grids, a power-law fit estimates the zero-spacing gain and convergence order.
The accompanying JSON retains the values and fit status. `--analyse-only STUDY`
can regenerate these products from the completed subset while other runs continue.

Install `tools/energy_analysis/requirements.txt` and the extra plotting package
`matplotlib`. Run from the repository root:

```sh
python examples/tokamak-2D-driftplane/2D-energy/publication_scan.py --help
python -m pytest examples/tokamak-2D-driftplane/2D-energy
```

This example is outside CTest. The legacy `2D-energy` regression remains in
`tests/integrated/2D-energy`; total ion energy conservation gates are the four
`ion-energy-*` integrated tests.

The input template is `publication-input.inp.template`. See the
[Sphinx comparison](https://hermes3.readthedocs.io/en/latest/examples.html#subsec-2d-energy-comparison)
for the energy definitions and reproduced figures.
