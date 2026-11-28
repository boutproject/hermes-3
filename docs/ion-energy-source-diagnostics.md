# Ion momentum and energy source diagnostics

Enable source attribution for a diagnostic window with:

```ini
[hermes]
diagnose_component_sources = true

[evolve_vorticity_ion_energy]
diagnose = true

[polarisation_drift]
diagnose = true

[d+]
diagnose = true

[e]
diagnose = true
```

The scheduler records each component's changes to ion `density_source`,
`momentum_source` and thermal `energy_source`. Differences are taken immediately
before and after that component's transform, so additive and replacing writes
are attributed correctly. Values present before the first component belong to
`initial_state`. No source values or time derivatives are changed.

For component C and ion S, output names are
`component_source_C_S_<quantity>`, with quantities:

| Quantity | Meaning | SI units after conversion |
| --- | --- | --- |
| `density` | Particle source increment, SN | m⁻³ s⁻¹ |
| `momentum` | Parallel momentum source increment, F | kg m⁻² s⁻² |
| `thermal_energy` | Thermal-energy source increment, Q | W m⁻³ |
| `momentum_work` | V F, using the completed velocity | W m⁻³ |
| `density_kinetic_work` | −½ A V² SN | W m⁻³ |

C is the component type when its instance has the same name, otherwise
`<type>_<instance>`. For example, friction uses `braginskii_friction` and an ion
momentum component uses `evolve_momentum_d+`. Sources from transform are recorded;
intrinsic time-derivative terms calculated in `finally()` have separate outputs.
The additional full-field storage and copying occur only when this option is
enabled. It applies to all declared ions; electron source attribution is not
included in this option.

The ion-energy component saves `Eion_source_C`, summed over ions, as

```text
Q_C + V F_C - 1/2 A V² SN_C
```

These are signed contributions to the energy rate. Conservative polarisation
already transports parallel kinetic energy in its energy flux. Its source
conversion is removed by the separate
`Eion_source_polarisation_conversion_correction` channel. Sum **all**
`Eion_source_*` channels, including that correction and `initial_state`, to
recover `Eion_src` to roundoff. The channels use final, boundary-conditioned
velocities, matching the energy source conversion used by the equation.

For each diagnosed momentum species S, `momentum_<term>_S` saves the signed
force density and `momentum_work_<term>_S` saves V times that force. Terms are
`exb`, `electric`, `parallel_advection`, `pressure`, `sources` and `other`.
Their sum is `ddt(NVS)`; `other` contains the remaining numerical, electromagnetic
and derivative-scaling contributions. `sources` is the aggregate external
momentum source; the scheduler outputs separate the component increments.

The full ion energy rate is

```text
ddt(Eion) = external_heating + Eion_src
            - DivEExB - DivEelectric - DivEparallel
            - DivEpol_work - DivEpol_transport
            - phi_Div_Je - electron_viscous_heating
```

`ion_external_heating`, `Eion_src`, the energy divergence fields,
`phi_Div_Je` and the optional viscous heating fields supply this budget.
`DivEpol_work = DivEpol_phi_work + DivEpol_psi_work`; do not add psi work a
second time. Save `ion_KE_parallel`, `ion_KE_perp`, `Eion`, density and temperature
to evaluate the recovered-energy deficit.

For a local parallel thermal-energy budget, recover the particle rate from the
saved area-integrated `pfd+_tot_ylow` flux. Its contribution at a physical cell is
`SN_parallel = -(flow[y+1] - flow[y]) / (J dx dy dz)`. The parallel kinetic rate is
`V momentum_parallel_advection_d+ - 1/2 A V² SN_parallel`; include
`momentum_work_pressure_d+` when comparing this with `-DivEparallel`.
Their difference is the parallel contribution to the recovered thermal-energy
rate. This comparison includes the numerical particle, momentum and energy
fluxes, which need not preserve positive thermal energy just because the total
energy flux is conservative. Additional transport and current work still enter
the complete thermal-energy budget.
For a shifted mesh, evaluate the face difference in field-aligned coordinates
and transform the rate back, or map the upper face into the cell's field-line
frame first. Subtracting unaligned neighbouring face values at the same stored
Z index is only valid for the Z integral or a Z-uniform state.

Comparing force work with total energy requires the particle kinetic conversion
as well as the energy fluxes. Pressure-force work is paired with enthalpy
transport. Electric work enters through the ion electrical flux and electron
current work, together with the perpendicular current/energy exchange.
An additional V F source for those terms would count work again. For friction
or particle exchange, the thermal source plus the two kinetic source-conversion
terms gives the component's contribution to total ion energy.

The output values are normalised and carry `conversion` metadata. Integrating a
rate over a local Z line gives a local diagnostic; rank-zero fields alone do not
provide a full-domain integral. The `balance_*` scalar ledgers are already MPI
reductions over the simulated domain.

For the downloaded target case, start from the complete `03` restart at
normalised time 1051 on the original 128-rank machine. The first saved deficit
at rank-zero `(x=2,y=4)` is bracketed by 1451 and 1551; the `04` restart at 1551
is already invalid there. Keep the baseline transport settings for this window,
enable the diagnostics above, and save at intervals of 0.1 near the onset. The
downloaded rank-zero restart alone cannot restart the complete domain.

At `(2,4)` and the corner `(2,2)`, compare each component's thermal source,
momentum work, density kinetic work and `Eion_source_C`, then the signed energy
divergences, particularly `DivEparallel`, `DivEpol_phi_work`,
`DivEpol_psi_work` and `DivEpol_transport`. Use volume weights for a Z-line
integral. The rate outputs describe the latest RHS evaluation; closely spaced
outputs are needed to compare their time integral with the evolved energy.
