"""Symbolic MMS identities and freshness of committed expression inputs."""

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


@pytest.mark.parametrize("variant", exact.VARIANTS)
def test_identities_and_periodicity(variant):
    report = exact.preflight(variant)
    assert report["positive_primitives"]
    assert report["elliptic_and_X_boundaries"] == "verified symbolically"
    for name in (*exact.EVOLVED, *exact.RECOVERED):
        for axis in (1, 2):
            points = [1.37, 2.91, 4.23, 0.73]
            points[axis] = 0
            first = exact.evaluate(name, *points, variant)
            points[axis] = 10
            assert np.allclose(
                first, exact.evaluate(name, *points, variant), atol=2e-13
            )


@pytest.mark.parametrize("variant", exact.VARIANTS)
def test_committed_inputs_are_current(tmp_path, variant):
    generated = tmp_path / "BOUT.inp"
    mms.generate(generated, variant)
    assert generated.read_bytes() == mms.input_path(variant).read_bytes()
    options = BoutOptionsFile(str(generated))
    assert options["solver"].get_bool("mms")
    for field in exact.EVOLVED:
        assert "source" in options[field] and "solution" in options[field]
    assert set(mms.REFERENCE_KEYS.values()) <= set(options["mms"].keys())


def test_common_temperature_single_ion_limits():
    ni, nj, psi = s.symbols("ni nj psi", positive=True)
    temperature = exact.common_temperature(psi, ni, nj)
    assert s.simplify(temperature.subs(nj, 0) - psi / ni) == 0
    assert s.simplify(temperature.subs(ni, 0) - 4 * psi / (3 * nj)) == 0
