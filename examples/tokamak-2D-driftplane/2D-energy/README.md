# Drift-plane energy scan

The scan compares matched legacy pressure and total ion energy formulations.
Both include `classical_diffusion`, with a constant normalised diffusivity
`custom_D=1.19` (17 m²/s for this setup), configurable with `--custom-D`.
It runs local cases and generates refinement tables and figures. See
[scan-study.md](scan-study.md) for the full study.
New studies are saved in `results/interchange-TIMESTAMP` beside the script;
`--output-root` selects another parent directory.

For weak diffusion and interchange dynamics, use `--custom-D 0.001`
(0.0143 m²/s); the historical default strongly damps the perturbation. The
initial pressure now has constant radial wall buffers and a smooth taper,
removing the pressure-model vorticity sheets caused by the previous profile.
This changes the initial state: start a new study to compare the formulations.

`ion-thermal-gain-<boundary>.pdf` plots the signed full-simulation gain
`U_i(T)-U_i(0)` against mesh spacing on a logarithmic axis. Available pressure
and ion-energy runs share the plot; with at least three suitable ion-energy
grids, a power-law fit estimates the zero-spacing gain and convergence order.
The accompanying JSON retains the values and fit status. `--analyse-only STUDY`
can regenerate these products from the completed subset while other runs continue.

Install `tools/energy_analysis/requirements.txt` and the extra plotting package
`matplotlib`. Run from the repository root:

```sh
python examples/tokamak-2D-driftplane/2D-energy/run_scan.py --help
python examples/tokamak-2D-driftplane/2D-energy/run_scan.py --smoke --custom-D 0.001
python -m pytest examples/tokamak-2D-driftplane/2D-energy
```

This example is outside CTest. `ion-energy-closed` is the small two-rank
conservation gate; `ion-energy-mms` is the optional 3D spatial-accuracy study.

The input template is `scan-input.inp.template`. See the
[Sphinx comparison](https://hermes3.readthedocs.io/en/latest/examples.html#subsec-2d-energy-comparison)
for the energy definitions and reproduced figures.
