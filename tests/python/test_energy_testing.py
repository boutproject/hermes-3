"""Acceptance thresholds must use resolved run tolerances and retain provenance."""

import json
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr
from boutdata.data import BoutOptions
from energy_analysis import AuditError
from energy_analysis.testing import Study, resolve_limit
from energy_analysis.walls import absolute_wall_flux


def run_with_tolerance(rtol):
    options = BoutOptions()
    options["accuracy"] = str(2 * rtol)
    options.getSection("solver")["rtol"] = "accuracy / 2"
    return SimpleNamespace(data=xr.Dataset(attrs={"options": options}))


def test_solver_relative_limit_uses_resolved_options():
    criterion = {"floor": 1e-10, "rtol_factor": 1000}
    assert resolve_limit(criterion, run_with_tolerance(1e-8)) == pytest.approx(1e-5)
    assert resolve_limit(criterion, run_with_tolerance(1e-14)) == 1e-10


@pytest.mark.parametrize("criterion", [0, 1e-11, {"floor": 1e-12}])
def test_roundoff_limit_cannot_resolve_below_floor(criterion):
    with pytest.raises(AuditError, match="below 1e-10"):
        resolve_limit(criterion)


def test_minimum_thresholds_are_exempt():
    assert resolve_limit(1e-12, minimum=True) == 1e-12
    assert resolve_limit(1.8, minimum=True) == 1.8


def test_mms_truncation_errors_are_exempt_from_roundoff_floor():
    assert resolve_limit(1e-14, truncation=True) == 1e-14


def test_solver_relative_limit_requires_a_run():
    with pytest.raises(AuditError, match="requires a run"):
        resolve_limit({"floor": 1e-10, "rtol_factor": 100})


@pytest.mark.parametrize("rtol", [0, -1e-8, float("nan")])
def test_invalid_run_tolerance_is_rejected(rtol):
    with pytest.raises((AuditError, ValueError)):
        resolve_limit({"floor": 1e-10, "rtol_factor": 100}, run_with_tolerance(rtol))


def test_effective_limit_is_written_with_check(tmp_path):
    criterion = {"floor": 1e-10, "rtol_factor": 1000}
    (tmp_path / "criteria.json").write_text(json.dumps({"restart": criterion}))
    args = SimpleNamespace(output_root=tmp_path)
    with Study(tmp_path, args) as study:
        study.check("restart", 1e-7, run=run_with_tolerance(1e-8))
    result = json.loads((study.output / "result.json").read_text())
    assert result["passed"]
    assert result["checks"][0]["limit"] == pytest.approx(1e-5)


def test_minimum_check_keeps_failure_visible(tmp_path):
    (tmp_path / "criteria.json").write_text(json.dumps({"flux": 1e-6}))
    with (
        pytest.raises(AuditError, match="Criteria failed"),
        Study(tmp_path, SimpleNamespace(output_root=tmp_path)) as study,
    ):
        study.check("flux", 0, minimum=True)
    result = json.loads((study.output / "result.json").read_text())
    assert not result["passed"]
    assert result["checks"][0]["comparison"] == ">="


def test_absolute_wall_flux_retains_cancelling_faces_and_excludes_y_guards():
    values = np.zeros((2, 8, 3, 2))
    values[:, [1, 5], 1, :] = 0.5
    values[:, :, [0, 2], :] = 100
    transport = xr.Dataset({"flux_x": (("t", "x", "y", "z"), values)})
    run = SimpleNamespace(mxg=2, myg=1)
    # Equal oriented fluxes enter one wall and leave the other, giving zero
    # net flux. The absolute inventory exchange remains finite.
    np.testing.assert_allclose(absolute_wall_flux(run, transport), 2)
