"""Failed MMS gates write a result before any subsequent solver launch."""

import json
import sys
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr
from energy_analysis.dataset import AuditError
from energy_analysis.testing import Study

RUNNER = (
    Path(__file__).resolve().parents[1] / "integrated/ion-energy-parallel-mms/runtest"
)


@pytest.mark.parametrize(
    "failure, expected",
    [
        ("physics", ["parallel-8"]),
        ("coarse_order", ["parallel-8", "parallel-16"]),
        ("field_error", ["parallel-8", "parallel-16"]),
        ("fine_order", ["parallel-8", "parallel-16", "parallel-32"]),
        ("temporal", ["parallel-8", "parallel-16", "parallel-32", "parallel-16-tight"]),
        (
            None,
            [
                "parallel-8",
                "parallel-16",
                "parallel-32",
                "parallel-16-tight",
                "parallel-8-decomposition",
            ],
        ),
    ],
)
def test_validation_stops_before_next_launch(tmp_path, monkeypatch, failure, expected):
    loader = SourceFileLoader("parallel_mms_workflow", str(RUNNER))
    runner = module_from_spec(spec_from_loader(loader.name, loader))
    loader.exec_module(runner)
    launched = []

    class SolverFreeStudy(Study):
        def run(self, label, overrides, **kwargs):
            launched.append(label)
            self.result["cases"][label] = {"outcome": "completed"}
            delta = 1e-6 if failure == "temporal" and label.endswith("tight") else 0
            field = xr.DataArray(np.full((2, 2), delta), dims=("t", "x"))
            run = SimpleNamespace(
                data=SimpleNamespace(metadata={"NXPE": 1, "NYPE": 1}, Pe=field + 1),
                field=field,
                interior=lambda value: value,
                integrate=lambda value: float(value.sum()),
            )
            self.current_run = run
            return run

    def physics(study, label, run):
        study.check(
            label + "/continuity",
            1 if failure == "physics" else 0,
            "continuity",
            truncation=True,
        )
        size = int(label.split("-")[1])
        error = 8e-6 * (8 / size) ** 2
        if failure == "coarse_order" and size == 16:
            error = 4e-6
        if failure == "field_error" and size == 16:
            error = 2e-5
        if failure == "field_error" and size == 8:
            error = 8e-5  # Passing order isolates the fine-grid field-error cap.
        if failure == "fine_order" and size == 32:
            error = 4e-6
        errors = {"NVe": {"l2": error, "linf": error}}
        study.result["cases"][label]["errors"] = errors
        return errors

    monkeypatch.setattr(runner, "Study", SolverFreeStudy)
    monkeypatch.setattr(runner, "check_physics", physics)
    monkeypatch.setattr(runner, "REFERENCES", {"NVe": "NVe:solution"})
    monkeypatch.setattr(runner, "comparison_fields", lambda run: {"NVe": run.field})
    monkeypatch.setattr(sys, "argv", [str(RUNNER), "--output-root", str(tmp_path)])
    if failure:
        with pytest.raises(AuditError, match="before further validation"):
            runner.main()
    else:
        runner.main()
    assert launched == expected
    [result_path] = tmp_path.glob("*/result.json")
    result = json.loads(result_path.read_text())
    assert result["passed"] is (failure is None)
    assert set(result["cases"]) == set(expected)
    assert any(not c["passed"] for c in result["checks"]) is (failure is not None)
    if failure == "field_error":
        assert [c["name"] for c in result["checks"] if not c["passed"]] == [
            "error/16/NVe"
        ]
