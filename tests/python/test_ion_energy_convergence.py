"""Nonlinear self-convergence of the fixed-wall publication drift plane."""

import json

import numpy as np
import pytest
from energy_analysis import open_run
from energy_analysis.inventory import energy_densities
from ion_energy_publication import run_case


@pytest.mark.slow
def test_publication_total_energy_self_convergence(hermes_executable, tmp_path):
    rows = []
    for n in (16, 32, 64):
        folder = run_case(
            hermes_executable, tmp_path / str(n), n, endpoint=5000, output_step=100
        )
        case = open_run(folder)
        try:
            ds = case.data
            assert float(ds.t[-1]) == 5000
            np.testing.assert_allclose(np.diff(ds.t), 100, atol=1e-12, rtol=0)
            kinetic = case.interior(energy_densities(case).ion_perpendicular.isel(t=-1))
            volume = case.interior(ds.J * ds.dx * ds.dy * ds.dz)
            energy = kinetic * volume
            total = float(energy.sum()) * case.energy_scale * case.length_scale**3
            # Exact cell/strip overlaps keep the physical outer 10% fixed.
            left = np.arange(n) / n
            overlap = np.maximum(0, left + 1 / n - np.maximum(left, 0.9)) * n
            radial = np.asarray(energy.sum(("y", "z")))
            share = float(np.dot(radial, overlap) / radial.sum())
            rows.append({"n": n, "kinetic_J": total, "outer_share": share})
        finally:
            case.data.close()
    (tmp_path / "convergence.json").write_text(json.dumps(rows, indent=2) + "\n")
    differences = abs(np.diff([r["kinetic_J"] for r in rows]))
    assert differences[0] > 1e-12
    assert differences[1] < 0.5 * differences[0], (rows, differences)
    shares = np.array([r["outer_share"] for r in rows])
    assert np.ptp(shares) < 0.005, rows
    assert shares[-1] <= shares[0] + 0.001, rows
