"""Symbolic validation of the single optional 3D slab MMS."""

import sys
from pathlib import Path

import numpy as np
import pytest

s = pytest.importorskip("sympy")
from boutdata.data import BoutOptionsFile

HERE = Path(__file__).resolve().parents[1] / "integrated/ion-energy-mms"
sys.path.insert(0, str(HERE))
import manufactured as exact
import mms


def test_continuum_identities_and_nonzero_terms():
    report = exact.preflight()
    assert report["positive_primitives"]
    assert report["g13"] == 0.25
    assert report["viscosity"] == 0.02
    assert report["magnetic_curvature"] == 0.1
    assert min(report["sample_term_magnitudes"].values()) > 1e-12


def test_periodic_y_and_z_for_evolved_and_recovered_fields():
    for name in (*exact.EVOLVED, *exact.RECOVERED):
        for axis in (1, 2):
            points = [1.37, 2.91, 4.23, 0.73]
            points[axis] = 0
            first = exact.evaluate(name, *points)
            points[axis] = float(exact.LENGTH)
            np.testing.assert_allclose(first, exact.evaluate(name, *points), atol=2e-13)


def test_generated_input_and_boundary_contract(tmp_path):
    path = tmp_path / "BOUT.inp"
    mms.generate(path)
    options = BoutOptionsFile(str(path))
    assert options["solver"].get_bool("mms")
    assert options["e:type"] == "quasineutral, evolve_pressure, evolve_momentum"
    assert options["Ve:bndry_xin"] == "dirichlet_o2(mms:electron_velocity)"
    for name in exact.EVOLVED:
        assert "source" in options[name] and "solution" in options[name]
    assert float(options["mesh"].evaluate_scalar("g13")) == 0.25
    assert (
        float(options["evolve_vorticity_ion_energy"].evaluate_scalar("viscosity"))
        == 0.02
    )
    assert float(options["mesh"].evaluate_scalar("bxcvz")) != 0


def test_sheared_slab_has_closed_viscous_faces_and_curved_pressure_walls():
    fields = exact.model()[0]
    H = exact.metric()
    assert H.det() == 1 and H[0, 0] > 0
    for face in (0, exact.LENGTH):
        assert exact.zero(fields["phi"].subs(exact.X, face))
        assert not exact.zero(s.diff(fields["phi"], exact.X).subs(exact.X, face))
        assert exact.zero(s.diff(fields["psi"], exact.X).subs(exact.X, face))
        assert not exact.zero(s.diff(fields["psi"], exact.X, 2).subs(exact.X, face))
        assert exact.zero(s.diff(fields["Vort"], exact.X).subs(exact.X, face))
        assert exact.zero((H * exact.grad(fields["Vort"]))[0].subs(exact.X, face))


def test_common_temperature_single_ion_limits():
    ni, nj, psi = s.symbols("ni nj psi", positive=True)
    temperature = exact.common_temperature(psi, ni, nj)
    assert s.simplify(temperature.subs(nj, 0) - psi / ni) == 0
    assert s.simplify(temperature.subs(ni, 0) - 4 * psi / (3 * nj)) == 0
