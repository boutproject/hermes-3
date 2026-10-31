# Two-ion MMS with a skew perpendicular metric

This focused variant of [the slab MMS](../3D-ion-energy-mms/README.md) uses
`g13=1/4`, `g33=17/16` and the consistent inverse metric, with `J=B=1`.
X remains nonperiodic; Y and Z remain periodic, without sheaths.

The exact electrostatic potential is zero. The pressure potential has a
radial factor `(1-cos(2*pi*X/10))^3`, which vanishes to sixth order at X
faces. This closes the mixed normal flux and satisfies the existing pressure
and auxiliary-potential boundaries. Vorticity and total energy are derived
from the continuum mixed operator and kinetic energy, rather than from
discrete operators. Parallel density, momentum and electron pressure
dynamics remain active. The baseline separately exercises Reynolds and ExB
transport.

Generate all inputs with `.venv/bin/python
tests/integrated/3D-ion-energy-mms/mms.py --all --check`, then run
`.venv/bin/python tests/integrated/3D-ion-energy-mms-skew/runtest` from the
repository root. The shared runner supports resolution, tolerance, MPI,
resume and analysis controls; reports include all evolved variables,
recovered pressures and potential, full-domain norms and X boundary errors.
