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
    for name in (*exact.evolved(variant), *exact.RECOVERED):
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
    for field in exact.evolved(variant):
        assert "source" in options[field] and "solution" in options[field]
    assert set(mms.REFERENCE_KEYS.values()) <= set(options["mms"].keys())
    if variant in exact.PARALLEL_VARIANTS:
        assert "evolve_momentum" in str(options["e:type"])
        assert "zero_current" not in str(options["e:type"])
        assert options["sound_speed:electron_dynamics"] == "true"


def test_parallel_variants_share_current_and_parallel_primitives():
    smooth, *_ = exact.model("parallel_current")
    wall, *_ = exact.model("parallel_current_wall")
    for name in ("Ni", "Nj", "NVi", "NVj", "NVe", "Te", "Pe", "Jpar", "DivJpar", "phi"):
        assert s.simplify(smooth[name] - wall[name]) == 0, name
    assert s.simplify(smooth["psi"] - wall["psi"]) != 0


@pytest.mark.parametrize("variant", exact.PARALLEL_VARIANTS)
def test_parallel_potential_curvature_at_both_walls(variant):
    fields, *_ = exact.model(variant)
    for face in (0, exact.LENGTH):
        for name in ("phi", "psi"):
            curvature = s.simplify(s.diff(fields[name], exact.X, 2).subs(exact.X, face))
            assert (curvature == 0) is (variant == "parallel_current" or name == "phi")
        for name in ("Ni", "Nj", "NVi", "NVj", "Tion", "Pe"):
            assert s.simplify(s.diff(fields[name], exact.Z).subs(exact.X, face)) == 0


@pytest.mark.parametrize("variant", exact.PARALLEL_VARIANTS)
def test_parallel_profiles_separate_reynolds_wall_error(variant):
    from types import SimpleNamespace

    import xarray as xr
    from energy_analysis.fluxes import fill_wall_guards
    from energy_analysis.transport import (
        face_divergence,
        reynolds_face_fluxes,
        xppm_faces,
    )
    from operator_convergence import expressions, sample

    fields, _, continuum = expressions(variant)
    errors = []
    for n in (16, 32, 64):
        h = float(exact.LENGTH) / n
        x, z = (np.arange(n + 4) - 1.5) * h, (np.arange(n) + 0.5) * h
        values = {
            name: sample(fields[name], x, z, 0.37, 0.4)
            for name in ("phi", "psi", "Vort")
        }
        wall = SimpleNamespace(mxg=2)
        phi = fill_wall_guards(values["phi"], wall, "zero")
        psi = fill_wall_guards(values["psi"], wall, "pressure")
        omega = fill_wall_guards(values["Vort"], wall, "neumann")
        geometry = xr.Dataset(
            {
                "J": 1.0,
                "Bxy": 1.0,
                "dx": h,
                "dy": 1.0,
                "dz": h,
                "g11": 1.0,
                "g13": 0.0,
                "g33": 1.0,
            }
        )
        faces = reynolds_face_fluxes(2.0, omega, phi, psi, geometry)
        advective = xppm_faces(omega, phi, geometry)
        remainder = {dim: faces[dim] - advective[dim] for dim in faces}
        error = np.asarray(
            face_divergence(remainder, geometry)
            - sample(continuum["R"], x, z, 0.37, 0.4)
        )[2:-2]
        interior = error[n // 4 : 3 * n // 4]
        errors.append(
            {
                "global_l2": np.sqrt(np.mean(error**2)),
                "interior_l2": np.sqrt(np.mean(interior**2)),
                "wall_max": np.max(abs(error[[0, -1]])),
            }
        )
    global_orders = np.log2(
        [
            errors[0]["global_l2"] / errors[1]["global_l2"],
            errors[1]["global_l2"] / errors[2]["global_l2"],
        ]
    )
    interior_orders = np.log2(
        [
            errors[0]["interior_l2"] / errors[1]["interior_l2"],
            errors[1]["interior_l2"] / errors[2]["interior_l2"],
        ]
    )
    assert np.all(interior_orders > 1.8), errors
    if variant == "parallel_current":
        assert np.all(global_orders > 1.8), errors
        assert errors[-1]["wall_max"] < errors[-2]["wall_max"] / 3
    else:
        assert global_orders[-1] < 1.0, errors
        assert errors[-1]["wall_max"] > 0.8 * errors[-2]["wall_max"], errors


def test_common_temperature_single_ion_limits():
    ni, nj, psi = s.symbols("ni nj psi", positive=True)
    temperature = exact.common_temperature(psi, ni, nj)
    assert s.simplify(temperature.subs(nj, 0) - psi / ni) == 0
    assert s.simplify(temperature.subs(ni, 0) - 4 * psi / (3 * nj)) == 0
