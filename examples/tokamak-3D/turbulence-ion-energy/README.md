# Ion total-energy turbulence

This example evolves ion density, ion total energy, electron pressure and both
parallel momenta. Electron pressure and ion external heating use the same
core mask. Ion heating is `1.5 * Pe:source` because `Pe` is a pressure, while
`external_heating` is an energy source in W/m³.

Install the `total-energy` branch of xhermes from this checkout, then launch
into a new directory and analyse the results separately:

```sh
python -m pip install -e external/xhermes
python examples/tokamak-3D/turbulence-ion-energy/run.py data --mz 16 --nproc 6
python examples/tokamak-3D/turbulence-ion-energy/analyse.py data
```

The launcher uses boutdata's `launch_safe`. Use `--grid` and `--exe` to select
other files. `--diagnose` enables the additional fields and saved RHS integrals;
otherwise content rates are obtained by differencing the outputs. A short
check can use `--mz 1 --nout 4 --timestep 0.02`.
`--set OPTION=VALUE` passes additional BOUT options. Use `--restart-from data
--nout 0 --diagnose` to evaluate a saved state without advancing it, keeping
the original MZ and number of ranks. `--input` selects another input file.

Analysis writes final-record tables, complete histories and PDF figures to
`data/balance`. It uses the actual saved source powers and particle rates,
including source masks and time dependence. The nominal `heating` and `flux`
settings depend on the mesh-specific `shape_factor`; they are not substituted
for measured domain integrals. Results cover the simulated toroidal sector
(`zperiod=5` here), including when MZ=1.

The energy table starts with `input_power`, `output_power`, `content_rate`,
and `power_imbalance = input_power - output_power - content_rate`. The imbalance
excludes bulk terms, so a nonzero value exposes a conservation defect even when
the equation ledger closes. It then lists ion/electron inputs, radial powers,
sheath powers and transport details, followed by the bulk work breakdown.
The figures show content, input/output/rate, sheath losses and residuals.
This grid has two sheath boundaries, at the beginning and end of Y. Four
target channels are supported for double-null grids.
`sheath_output` sums the targets, including their signed ion/electron electrical
power; `sheath_electric_*` shows this part of each target's total. Older dumps
may lack the electrical target split. `remaining_output` is the remaining
transport integral, including `out_electric_perpendicular` and any volume
defects in transport. Inner and outer radial polarisation powers are listed separately.
`bulk` retains the net non-transport equation terms, including source
conversions, mechanical work and discrete conservation defects.
It must not be interpreted as external
heating. `boundary_residual` shows input minus transport losses minus content
rate; `residual` also includes the bulk terms. Incomplete older outputs retain
their inventories and known inputs but leave unavailable balances as NaN.

Parallel numerical viscous heating is included for both species: ion heating
is retained by total-energy recovery, and `e:numerical_viscous_heating=true`
returns the electron correction to pressure. To compare with electron heating
disabled, use `--set e:numerical_viscous_heating=false`. The correction includes
sheath momentum-flux terms and can be negative. Perpendicular
ion viscosity and collisional frictional heating are separate mechanisms.
The Braginskii parallel ion/electron viscosity components are not enabled.
The `bulk_*` rows show the terms summed into `bulk`; diagnostics also show
the available parallel numerical-heating correction when it is disabled.
The electrical mismatch rows compare the actual work terms with the reported
electrical transport. Electron parallel output uses the cross-product sheath
flux implied by the centred force/current operators. Its difference from the
average of cell/guard products is a boundary accounting correction; the
evolution equations are unchanged. Older dumps use their original accounting.
