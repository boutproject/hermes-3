"""MMS gate sequencing, selected norm regions and publication failure reports."""

import json
import subprocess
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

RUNNER = Path(__file__).resolve().parents[1] / "integrated/ion-energy-mms/runtest"


def test_runner_does_not_require_symbolic_generator():
    code = """
import importlib.abc
import runpy
import sys
class NoSympy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, *args):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise ImportError('SymPy disabled for this runtime check')
sys.meta_path.insert(0, NoSympy())
sys.argv = [sys.argv[1], '--help']
runpy.run_path(sys.argv[0], run_name='__main__')
"""
    result = subprocess.run(
        [sys.executable, "-c", code, str(RUNNER)],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "--levels" in result.stdout and "--variant" not in result.stdout


@pytest.fixture
def runner():
    loader = SourceFileLoader("parallel_mms_workflow", str(RUNNER))
    module = module_from_spec(spec_from_loader(loader.name, loader))
    loader.exec_module(module)
    return module


def fake_workflow(runner, monkeypatch, failure=None):
    launched, plots = [], []

    class SolverFreeStudy(Study):
        def run(self, label, overrides, **kwargs):
            launched.append(label)
            self.result["cases"][label] = {"outcome": "completed"}
            if failure == "execution" and label.endswith("-32"):
                self.result["cases"][label]["outcome"] = "stalled"
                raise AuditError("Run stalled")
            delta = (
                1e-6
                if (
                    failure == "temporal"
                    and label.endswith("tight")
                    or failure == "decomposition"
                    and label.endswith("decomposition")
                )
                else 0
            )
            field = xr.DataArray(np.full((2, 2), delta), dims=("t", "x"))
            run = SimpleNamespace(
                data=SimpleNamespace(metadata={"NXPE": 1, "NYPE": 1}, Pe=field + 1),
                field=field,
                interior=lambda value: value,
                integrate=lambda value: float(value.sum()),
            )
            self.current_run = run
            return run

    def physics(study, label, run, variant):
        size = int(label.split("-")[1])
        failed = failure == "physics" or failure == "late_physics" and size == 32
        study.check(label + "/continuity", int(failed), "continuity", truncation=True)
        error = 8e-6 * (8 / size) ** 2
        if failure == "coarse_order" and size == 16:
            error = 4e-6
        if failure == "large_converging_error":
            error *= 1e6
        if failure == "fine_order" and size == 32:
            error = 4e-6
        chi = 0.016 * (8 / size) ** 2
        if failure == "chi_nonconverging" and size == 16:
            chi = 0.016
        if failure == "chi_nonfinite" and size == 16:
            raise AuditError("Nonfinite chi before further validation")
        errors = {
            name: dict.fromkeys(("l2", "linf", "l2_interior", "linf_interior"), value)
            for name, value in (("NVe", error), ("ion_chi", chi))
        }
        study.result["cases"][label]["errors"] = errors
        return errors

    def plot(output, variant, errors):
        plots.append(list(errors))
        (output / "convergence.csv").write_text("mock publication data\n")
        return ["convergence.csv"]

    import plotting

    # CLI/workflow tests do not depend on the optional plotting installation.
    monkeypatch.setitem(sys.modules, "matplotlib", SimpleNamespace(use=lambda _: None))
    monkeypatch.setattr(plotting, "plot_convergence", plot)
    monkeypatch.setattr(runner, "Study", SolverFreeStudy)
    monkeypatch.setattr(runner, "check_physics", physics)
    monkeypatch.setattr(
        runner,
        "references",
        lambda variant: {"NVe": "NVe:solution", "ion_chi": "mms:current_potential"},
    )
    monkeypatch.setattr(
        runner,
        "comparison_fields",
        lambda run, variant: {"NVe": run.field, "ion_chi": run.field},
    )
    return launched, plots


@pytest.mark.parametrize(
    "failure, suffixes",
    [
        ("physics", ["8"]),
        ("coarse_order", ["8", "16"]),
        ("chi_nonconverging", ["8", "16"]),
        ("chi_nonfinite", ["8", "16"]),
        ("fine_order", ["8", "16", "32"]),
        ("temporal", ["8", "16", "32", "16-tight"]),
        ("decomposition", ["8", "16", "32", "16-tight", "8-decomposition"]),
        (None, ["8", "16", "32", "16-tight", "8-decomposition"]),
    ],
)
def test_validation_stops_before_next_launch(
    tmp_path, monkeypatch, runner, failure, suffixes
):
    launched, plots = fake_workflow(runner, monkeypatch, failure)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(RUNNER),
            "--output-root",
            str(tmp_path),
            "--levels",
            "3",
            "--verify",
        ],
    )
    if failure:
        with pytest.raises(AuditError):
            runner.main()
    else:
        runner.main()
    expected = ["parallel_current_wall-" + suffix for suffix in suffixes]
    assert launched == expected
    assert plots == []
    [result_path] = tmp_path.glob("*/result.json")
    result = json.loads(result_path.read_text())
    assert result["passed"] is (failure is None)
    assert set(result["cases"]) == set(expected)
    assert any(not c["passed"] for c in result["checks"]) is (
        failure is not None and failure != "chi_nonfinite"
    )
    assert not any(c["name"].startswith("error/") for c in result["checks"])


def test_two_levels_allow_large_converging_chi(tmp_path, monkeypatch, runner):
    launched, _ = fake_workflow(runner, monkeypatch)
    monkeypatch.setattr(
        sys, "argv", [str(RUNNER), "--output-root", str(tmp_path), "--verify"]
    )
    runner.main()
    assert launched == [
        "parallel_current_wall-8",
        "parallel_current_wall-16",
        "parallel_current_wall-16-tight",
        "parallel_current_wall-8-decomposition",
    ]
    [path] = tmp_path.glob("*/result.json")
    result = json.loads(path.read_text())
    assert result["passed"]
    assert (
        result["cases"]["parallel_current_wall-16"]["errors"]["ion_chi"]["l2"] == 0.004
    )


@pytest.mark.parametrize(
    "options, levels",
    [
        ([], 2),
        (["--plot"], 4),
        (["--plot", "--levels", "5"], 5),
        (["--levels", "3"], 3),
    ],
)
def test_level_defaults(monkeypatch, runner, options, levels):
    monkeypatch.setitem(sys.modules, "matplotlib", SimpleNamespace(use=lambda _: None))
    monkeypatch.setattr(sys, "argv", [str(RUNNER), *options])
    assert runner.parse_args().levels == levels


def test_invalid_levels_rejected(monkeypatch, runner):
    monkeypatch.setattr(sys, "argv", [str(RUNNER), "--levels", "1"])
    with pytest.raises(SystemExit):
        runner.parse_args()


@pytest.mark.parametrize(
    "failure, sizes",
    [
        ("coarse_order", [8, 16, 32, 64]),
        ("late_physics", [8, 16, 32]),
        ("execution", [8, 16]),
    ],
)
def test_plot_failure_saves_completed_levels(
    tmp_path, monkeypatch, runner, failure, sizes
):
    launched, plots = fake_workflow(runner, monkeypatch, failure)
    monkeypatch.setattr(
        sys, "argv", [str(RUNNER), "--output-root", str(tmp_path), "--plot", "--verify"]
    )
    with pytest.raises(AuditError):
        runner.main()
    assert plots == [sizes]
    expected_launches = [8, 16, 32] if failure == "execution" else sizes
    assert launched == [f"parallel_current_wall-{size}" for size in expected_launches]
    [path] = tmp_path.glob("*/result.json")
    result = json.loads(path.read_text())
    assert not result["passed"]
    assert result["plots"] == ["convergence.csv"]
    assert (path.parent / "convergence.csv").is_file()
    if failure == "coarse_order":
        assert result["verification"] == "not run: spatial accuracy failed"
        assert set(result["orders"]) == {"16", "32", "64"}


def test_curved_wall_gates_interior_and_reports_global(tmp_path, runner):
    args = SimpleNamespace(output_root=tmp_path)
    coarse = {
        "ion_chi": {"l2": 0.02, "linf": 0.1, "l2_interior": 0.02, "linf_interior": 0.1}
    }
    fine = {
        "ion_chi": {
            "l2": 0.018,
            "linf": 0.1,
            "l2_interior": 0.005,
            "linf_interior": 0.025,
        }
    }
    runner.references = lambda variant: {"ion_chi": "mms:current_potential"}
    study = Study(runner.HERE, args)
    runner.check_refinement(study, coarse, fine, 16, "parallel_current_wall")
    assert all(c["passed"] for c in study.result["checks"])
    assert study.result["orders"]["16"]["ion_chi"]["linf"] == 0
    assert study.result["orders"]["16"]["ion_chi"]["linf_interior"] == 2


def test_mc_aware_parallel_limits(tmp_path, runner):
    study = Study(runner.HERE, SimpleNamespace(output_root=tmp_path))
    runner.references = lambda variant: {"ion_chi": "mms:current_potential"}
    coarse = {
        "ion_chi": dict.fromkeys(("l2", "linf", "l2_interior", "linf_interior"), 0.01)
    }
    fine = {
        "ion_chi": {
            "l2": 0.01 / 2**1.5,
            "linf": 0.005,
            "l2_interior": 0.01 / 2**1.5,
            "linf_interior": 0.005,
        }
    }
    runner.check_refinement(study, coarse, fine, 16, "parallel_current_wall")
    assert all(c["passed"] for c in study.result["checks"])
    assert [c["limit"] for c in study.result["checks"]] == [1.4, 0.9]


def test_completed_rate_subtracts_one_continuum_mms_source(monkeypatch, runner):
    field = xr.DataArray(np.ones((2, 3)), dims=("t", "x"), coords={"t": [0.0, 0.02]})
    ds = xr.Dataset(
        {
            "ddt(Vort)": 10 * field,
            "ion_current_rhs": 2 * field,
            "DivJviscous": 3 * field,
            "ion_phi_dissipation": 4 * field,
        }
    )
    calls = []

    def forcing(run, key, time):
        calls.append((key, time))
        return field.isel(t=0, drop=True)

    monkeypatch.setattr(runner, "evaluate", forcing)
    np.testing.assert_array_equal(
        runner.vorticity_rate_residual(SimpleNamespace(data=ds)), 0
    )
    assert calls == [("Vort:source", 0.0), ("Vort:source", 0.02)]


def test_magnetic_current_uses_slab_metric_and_physical_normalisation(
    monkeypatch, runner
):
    z = 2 * np.pi * np.arange(16) / 16
    sine = xr.DataArray(np.sin(z), dims="z")
    cosine = xr.DataArray(np.cos(z), dims="z")
    ds = xr.Dataset({"Pe": 2 + sine, "Bxy": 2.0, "dz": 2 * np.pi / 16})
    ds.attrs["options"] = {
        "mesh": SimpleNamespace(evaluate_scalar=lambda _: 0.0125),
        "hermes": SimpleNamespace(evaluate_scalar=lambda _: 1.0),
    }
    run = SimpleNamespace(data=ds, length_scale=2.0)
    monkeypatch.setattr(
        runner,
        "energy_densities",
        lambda _: SimpleNamespace(ion_thermal=1.5 * (1 + 0.2 * cosine)),
    )
    np.testing.assert_allclose(
        runner.diamagnetic_divergence(run),
        0.05 * (cosine - 0.2 * sine),
        atol=5e-16,
        rtol=0,
    )


def test_spatial_acceptance_uses_orders_without_absolute_error_caps(
    tmp_path, monkeypatch, runner
):
    launched, _ = fake_workflow(runner, monkeypatch, "large_converging_error")
    monkeypatch.setattr(
        sys, "argv", [str(RUNNER), "--output-root", str(tmp_path), "--verify"]
    )
    runner.main()
    [result_path] = tmp_path.glob("*/result.json")
    result = json.loads(result_path.read_text())
    assert result["passed"]
    assert result["cases"]["parallel_current_wall-16"]["errors"]["NVe"]["l2"] == 2
    assert not any(c["name"].startswith("error/") for c in result["checks"])
    assert launched[-1].endswith("decomposition")
