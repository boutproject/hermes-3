# Drift-plane energy publication scan

The scan compares matched legacy pressure and total ion energy formulations.
It generates refinement tables and figures and supports preparing and importing
remote runs. See [publication-study.md](publication-study.md) for the full study.

Install `tools/energy_analysis/requirements.txt` and the extra plotting package
`matplotlib`. Run from the repository root:

```sh
python examples/tokamak-2D-driftplane/2D-energy/publication_scan.py --help
python -m pytest examples/tokamak-2D-driftplane/2D-energy
```

This example is outside CTest. The legacy `2D-energy` regression remains in
`tests/integrated/2D-energy`; total ion energy conservation gates are the four
`ion-energy-*` integrated tests.
