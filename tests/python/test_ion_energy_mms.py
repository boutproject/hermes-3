"""Contracts for analytic MMS sources, boundaries and convergence failures."""

import importlib.machinery
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
import sympy as s
from boutdata.data import BoutOptionsFile
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parents[1] / "integrated/3D-ion-energy-mms"
sys.path.insert(0, str(HERE))
import analysis as mms_analysis
import manufactured as exact
import mms

loader = importlib.machinery.SourceFileLoader("ion_mms_runner", str(HERE / "runtest"))
spec = importlib.util.spec_from_loader(loader.name, loader)
runner = importlib.util.module_from_spec(spec)
loader.exec_module(runner)


@pytest.mark.parametrize("variant", exact.VARIANTS)
def test_analytic_current_and_physical_boundaries(variant):
    report = exact.preflight(variant)
    assert report["positive_primitives"]
    assert report["elliptic_and_X_boundaries"] == "verified symbolically"
    for name in (*exact.EVOLVED, *exact.RECOVERED):
        # Periodic Y/Z seam checks include all independently derived fields.
        for axis in (1, 2):
            points = [1.37, 2.91, 4.23, 0.73]
            points[axis] = 0
            first = exact.evaluate(name, *points, variant)
            points[axis] = 10
            assert np.allclose(
                first, exact.evaluate(name, *points, variant), atol=2e-13
            )


def test_common_temperature_single_ion_limits():
    ni, nj, psi = s.symbols("ni nj psi", positive=True)
    temperature = exact.common_temperature(psi, ni, nj)
    assert s.simplify(temperature.subs(nj, 0) - psi / ni) == 0
    assert s.simplify(temperature.subs(ni, 0) - 4 * psi / (3 * nj)) == 0


@pytest.mark.parametrize("variant", exact.VARIANTS)
def test_generated_options_reserve_sources_and_define_references(variant):
    options = BoutOptionsFile(str(mms.input_path(variant)))
    assert options["solver"].get_bool("mms")
    assert not options.get_bool("periodicX")
    assert int(options["mesh"].evaluate_scalar("jyseps1_1")) == -1
    assert options["mms:variant"] == variant
    for field in exact.EVOLVED:
        assert "source" in options[field] and "solution" in options[field]
        assert "**" not in options[field]["source"]
    assert len(set(mms.REFERENCE_KEYS.values())) == len(exact.RECOVERED)
    # A case-insensitive options reader must preserve both total and species p_i.
    assert set(mms.REFERENCE_KEYS.values()) <= set(options["mms"].keys())
    if variant == "skew":
        assert options["mesh"].evaluate_scalar("g13") == pytest.approx(0.25)
        H = exact.metric(variant)
        assert H.det() == 1


def test_output_completion_checks_every_rank_and_field(tmp_path):
    def write_rank(rank):
        with Dataset(tmp_path / f"BOUT.dmp.{rank}.nc", "w") as ds:
            ds.createDimension("t", 2)
            ds.createVariable("t", "f8", ("t",))[:] = [0, 0.1]
            for name in (
                *exact.EVOLVED,
                "phi",
                "Tion",
                *("E_" + v for v in exact.EVOLVED),
            ):
                ds.createVariable(name, "f8", ("t",))[:] = [1, 1]

    write_rank(0)
    assert runner.output_complete(tmp_path, 1, 0.1, 1)
    assert not runner.output_complete(tmp_path, 2, 0.1, 1)
    write_rank(1)
    assert runner.output_complete(tmp_path, 2, 0.1, 1)
    with Dataset(tmp_path / "BOUT.dmp.1.nc", "a") as ds:
        ds["phi"][-1] = np.nan
    assert not runner.output_complete(tmp_path, 2, 0.1, 1)
    with Dataset(tmp_path / "BOUT.dmp.1.nc", "a") as ds:
        ds["phi"][-1] = 1
        ds["t"][-1] = 0.05
    assert not runner.output_complete(tmp_path, 2, 0.1, 1)


def test_full_domain_norms_retain_boundary_cells_and_reject_nan():
    error = np.array([4.0, 0, 0, 0])
    volume = np.array([1.0, 2, 3, 4])
    assert mms_analysis.norms(error, volume) == pytest.approx((np.sqrt(1.6), 4))
    with pytest.raises(ValueError, match="Nonfinite"):
        mms_analysis.norms(np.array([np.nan]), np.ones(1))


def test_two_resolutions_cannot_silently_pass_criteria(tmp_path):
    with pytest.raises(ValueError, match="at least three resolutions"):
        mms_analysis.write_scan(
            tmp_path, [{}, {}], criteria={"default": {"l2_min_order": 1.5}}
        )


def test_nonconvergence_fails_even_with_three_completed_cases(tmp_path):
    reports = []
    for size, error in ((8, 1.0), (16, 0.5), (32, 0.4)):
        reports.append(
            {
                "errors": [
                    {
                        "case": str(size),
                        "field": "phi",
                        "time": 0.1,
                        "nx": size,
                        "ny": size,
                        "nz": size,
                        "l2": error,
                        "linf": error,
                    }
                ]
            }
        )
    failures = mms_analysis.write_scan(
        tmp_path, reports, criteria={"phi": {"l2_min_order": 1.5}}
    )
    assert len(failures) == 1 and "phi" in failures[0]
