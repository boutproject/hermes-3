# Radial ion energy walls

An orthogonal slab has physical radial walls and periodic Z, with hydrogen and
quasineutral electrons and no parallel flow. Six small cases cover `fixed`,
`free_slip` and `no_slip`, each with Neumann and Dirichlet AC conditions. They
check wall powers, particles, energy, vorticity and stored wall potentials from
ordinary outputs. An open AC wall without chi reports incomplete coverage.
The no-slip cases also test rigid-rotation decay and exact uniform-vorticity
drift; a 32-cell case measures second-order asymptotic amplitude convergence.
Fixed and no-slip runs have advancing restart round trips.
Their checkpoint segments also check diagnostics on/off against the matching
records of uninterrupted runs; every fresh run saves its initial state.

```sh
python tests/integrated/ion-energy-walls/runtest --exe build/hermes-3
```

The default uses 11 executions: seven budget/analytic references and four
checkpoint/restart segments. Separate restart segments are necessary for a
comparison with uninterrupted evolution. All results and limits appear in
`result.json` and `criteria.json`. The validated profile took 21 seconds.
No sheaths or parallel current are present.
