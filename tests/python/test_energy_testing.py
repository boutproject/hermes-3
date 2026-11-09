"""Acceptance thresholds must use resolved run tolerances and retain provenance."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr
from boutdata.data import BoutOptions
from energy_analysis import AuditError, testing
from energy_analysis.testing import Study, resolve_limit, run_progress
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


def test_progress_ignores_options_and_reads_saved_output_rows():
    log = "\tOption density_floor = 1e-7 (default)\n"
    log += "\tOption test = 1 (Minimum density floor)\n"
    log += "0.000e+00    1    0.01\n| Step 0\r2.000e+00    10    0.02\n"
    assert run_progress(log) == (2.0, None)
    assert run_progress(log + "Ion energy below floor\n") == (
        2.0,
        "Ion energy below floor",
    )


@pytest.mark.parametrize(
    ("events", "exit_time", "returncode", "stall", "timeout", "outcome", "last_time"),
    [
        ([], 1, 0, 120, 180, "completed", 0.0),
        ([], 1, 3, 120, 180, "failed", 0.0),
        ([(1, "Ion energy below floor\n")], 200, 0, 120, 180, "floored", 0.0),
        ([(1, "Density floor reached\n")], 1, 0, 120, 180, "floored", 0.0),
        ([], 200, 0, 120, 180, "stalled", 0.0),
        ([], 200, 0, 120, 30, "failed", 0.0),
        ([(110, "2.000e+00 10 0.1\n")], 250, 0, 120, 300, "stalled", 2.0),
        ([(110, "0.000e+00 10 0.1\n")], 200, 0, 120, 180, "stalled", 0.0),
    ],
)
def test_launch_records_monitored_outcomes(
    tmp_path,
    monkeypatch,
    events,
    exit_time,
    returncode,
    stall,
    timeout,
    outcome,
    last_time,
):
    # A controlled clock and output schedule test the actual 120-second stall
    # rule without depending on process startup or scheduler latency.
    clock = [0.0]
    monkeypatch.setattr(
        testing, "revision", lambda root: {"revision": "test", "dirty": False}
    )
    monkeypatch.setattr(testing.time, "monotonic", lambda: clock[0])

    def advance(seconds):
        clock[0] += seconds

    monkeypatch.setattr(testing.time, "sleep", advance)

    class Process:
        pid = 12345
        returncode = None

        def __init__(self, command, *, stdout, stderr, start_new_session):
            assert start_new_session
            self.log = stdout
            self.events = list(events)
            self.log.write("Option density_floor = 1e-7 (default)\n0.000e+00 1 0.01\n")
            self.log.flush()
            (Path(command[2]) / "BOUT.dmp.0.nc").touch()

        def poll(self):
            while self.events and self.events[0][0] <= clock[0]:
                self.log.write(self.events.pop(0)[1])
                self.log.flush()
            if clock[0] >= exit_time:
                self.returncode = returncode
            return self.returncode

        def wait(self, timeout=None):
            self.returncode = -15

    monkeypatch.setattr(testing.subprocess, "Popen", Process)
    signals = []
    monkeypatch.setattr(
        testing.os, "killpg", lambda pid, sig: signals.append((pid, sig))
    )
    (tmp_path / "criteria.json").write_text(json.dumps({"stall_seconds": stall}))
    executable = tmp_path / "hermes"
    executable.write_text("fake executable for provenance")
    args = SimpleNamespace(
        output_root=tmp_path, exe=executable, nproc=1, timeout=timeout
    )
    study = Study(tmp_path, args)
    if outcome == "completed":
        with study:
            folder = study.output / "case"
            folder.mkdir()
            (folder / "BOUT.inp").write_text("[solver]\nrtol=1e-7\n")
            study.launch(folder)
    else:
        with pytest.raises(AuditError), study:
            folder = study.output / "case"
            folder.mkdir()
            (folder / "BOUT.inp").write_text("[solver]\nrtol=1e-7\n")
            study.launch(folder)
    result = json.loads((study.output / "result.json").read_text())
    status = result["cases"]["case"]
    assert status["outcome"] == outcome
    assert status["last_output_time"] == last_time
    if outcome == "stalled":
        expected = 110 + stall if last_time == 2 else stall
        assert expected <= status["elapsed_seconds"] <= expected + 1
    if outcome == "stalled" or status.get("timed_out"):
        assert signals
    assert result["passed"] == (outcome == "completed")
    if outcome != "completed":
        assert not result["checks"][-1]["passed"]
    if timeout < stall:
        assert status["timed_out"]
