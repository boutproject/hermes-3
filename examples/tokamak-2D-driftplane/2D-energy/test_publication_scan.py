"""Independent analytic checks for publication estimators and field mapping."""

import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))

import numpy as np
import pytest
import xarray as xr
from energy_analysis import AuditError
from energy_analysis.convergence import (
    error_norms,
    global_differences,
    overlap_matrix,
    restrict,
)
from energy_analysis.pressure import legacy_operator


def grid(nx=20, nz=32):
    x, z = (np.arange(nx + 4) - 1.5) / nx, np.arange(nz) * 2 * np.pi / nz
    f = xr.DataArray(
        np.zeros((nx + 4, 1, nz)),
        dims=("x", "y", "z"),
        coords={"x": x, "y": [0], "z": z},
    )
    g = xr.Dataset(
        {
            n: xr.ones_like(f) * v
            for n, v in {
                "J": 1,
                "dx": 1 / nx,
                "dz": 2 * np.pi / nz,
                "g11": 1,
                "g33": 1,
                "g13": 0,
                "Bxy": 1,
            }.items()
        }
    )
    return f, g


def test_native_spectral_energy_and_radial_boundary_work():
    f, g = grid()
    f[:] = 0.3 * np.cos(3 * f.z) + 0.2 * f.x
    result = legacy_operator(2, f, g)
    interior = lambda v: v.isel(x=slice(2, -2))
    volume = g.J * g.dx * g.dz
    energy = float((interior(result.kperp * volume)).sum())
    # Linear radial field and an exactly represented Fourier mode.
    expected = 2 * np.pi * (0.5 * 2 * 0.2**2 + 0.25 * 2 * 0.3**2 * 3**2)
    assert energy == pytest.approx(expected, rel=1e-13)
    # -f L f/2 alone misses radial midpoint work.
    bulk = float(interior(-0.5 * f * result.divergence * volume).sum())
    face = lambda c, h: (
        (f.isel(x=c) + f.isel(x=h)) / 2 * 2 * (f.isel(x=h) - f.isel(x=c)) / (1 / 20)
    )
    boundary = float(((face(-3, -2) + face(2, 1)) * (2 * np.pi / 32)).sum())
    assert energy == pytest.approx(bulk + boundary / 2, abs=1e-13)


def test_native_nyquist_energy_is_retained():
    f, g = grid(nz=16)
    f[:] = np.cos(8 * f.z)
    result = legacy_operator(1, f, g)
    assert float(result.kperp.isel(x=5).mean()) == pytest.approx(32)


def test_periodic_shifted_faces_preserve_constant_and_integral():
    fine = np.arange(32) / 32
    coarse = np.arange(8) / 8
    r = overlap_matrix(fine, coarse, source_width=1 / 32, target_width=1 / 8, period=1)
    assert np.allclose(r @ np.ones(32), 1)
    values = np.sin(2 * np.pi * fine) + 2
    assert np.mean(r @ values) == pytest.approx(np.mean(values))
    # First coarse cell straddles the periodic seam.
    assert r[0, -1] > 0
    assert np.argmax(r[0]) != 1


def test_point_mapping_is_accurate_for_polynomial_and_fourier():
    sx = (np.arange(64) + 0.5) / 64
    sz = np.arange(64) / 64
    tx = (np.arange(8) + 0.5) / 8
    tz = np.arange(8) / 8
    field = sx[:, None] ** 3 * np.sin(6 * np.pi * sz)[None, :]
    mapped = restrict(field, sx, sz, tx, tz, lx=1, lz=1, representation="point")
    assert (
        np.max(abs(mapped - tx[:, None] ** 3 * np.sin(6 * np.pi * tz)[None, :])) < 1e-14
    )


def test_cell_average_mapping_preserves_constant_and_integral():
    sx = (np.arange(32) + 0.5) / 32
    sz = np.arange(32) / 32
    tx = (np.arange(8) + 0.5) / 8
    tz = np.arange(8) / 8
    field = np.ones((32, 32)) * 7
    mapped = restrict(field, sx, sz, tx, tz, lx=1, lz=1)
    assert np.allclose(mapped, 7)
    field = np.arange(32 * 32).reshape(32, 32)
    mapped = restrict(field, sx, sz, tx, tz, lx=1, lz=1)
    assert np.mean(mapped) == pytest.approx(np.mean(field))


def test_weighted_norms_keep_mean_errors_and_zero_references():
    result = error_norms(np.ones((2, 2)), np.zeros((2, 2)), np.ones((2, 2)))
    assert result == {"L2": 1, "Linf": 1, "relative_L2": None, "relative_Linf": None}
    result = error_norms(np.ones((2, 2)), np.ones((2, 2)), np.ones((2, 2)))
    assert result["L2"] == result["Linf"] == 0


def test_adjacent_orders_and_oscillatory_roundoff_sequences():
    rows = global_differences([1 + 4**-i for i in range(4)], 2)
    assert rows[0]["observed_order"] == pytest.approx(2)
    assert rows[1]["observed_order"] == pytest.approx(2)
    assert rows[-1]["finest_difference_J"] is None
    assert global_differences([1, 2, 1], 2)[0]["observed_order"] is None
    assert global_differences([1, 1, 1], 2)[0]["order_status"] == "roundoff dominated"
    assert global_differences([1, 2], 2)[0]["observed_order"] is None


def test_mismatched_boxes_fail():
    with pytest.raises(AuditError):
        overlap_matrix([0.25, 0.75], [0.4], source_width=0.5, target_width=1)


@pytest.fixture
def runner():
    source = Path(__file__).resolve().parent / "publication_scan.py"
    spec = importlib.util.spec_from_file_location("publication_scan", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_resolutions_guards_and_cadence(runner):
    args = runner.parser().parse_args(["--end-time", "120", "--outputs", "20"])
    config = runner.configuration(args)
    assert config["snapshot_times"] == [0, 30, 60, 120]
    assert config["outputs"] == 20
    text = runner.input_text(
        config, {"nx": 64, "nz": 64, "formulation": "total_energy", "boundary": "fixed"}
    )
    assert "nx = 68" in text
    assert "type = evolve_density\n" in text
    assert "[Ph+]" not in text
    assert "inner_boundary_flags" not in text
    assert "outer_boundary_flags" not in text
    pressure_text = runner.input_text(
        config, {"nx": 64, "nz": 64, "formulation": "pressure", "boundary": "fixed"}
    )
    assert "inner_boundary_flags = 16" in pressure_text
    assert "outer_boundary_flags = 16" in pressure_text
    args = runner.parser().parse_args(["--end-time", "120", "--outputs", "10"])
    assert runner.configuration(args)["outputs"] == 20


@pytest.mark.parametrize("coefficient", [1.19, 0.4])
def test_matched_inputs_include_constant_classical_diffusion(runner, coefficient):
    config = runner.configuration(
        runner.parser().parse_args(["--custom-D", str(coefficient), "--diagnose"])
    )
    for form in ("pressure", "total_energy"):
        text = runner.input_text(
            config, {"nx": 16, "nz": 16, "formulation": form, "boundary": "fixed"}
        )
        components = next(
            line for line in text.splitlines() if line.startswith("components =")
        )
        assert "classical_diffusion" in components
        assert (
            f"[classical_diffusion]\ncustom_D = {coefficient}\ndiagnose = true" in text
        )


@pytest.mark.parametrize("coefficient", ["0", "-1", "nan", "inf"])
def test_invalid_custom_diffusion_rejected(runner, coefficient):
    with pytest.raises(AuditError, match="custom_D must be finite and positive"):
        runner.configuration(runner.parser().parse_args(["--custom-D=" + coefficient]))


@pytest.mark.parametrize("order", [0.5, 2, 4])
@pytest.mark.parametrize("limit,amplitude", [(3e-7, 2e-6), (-3e-7, -2e-6)])
def test_thermal_gain_fit_recovers_signed_nonzero_limit(
    runner, order, limit, amplitude
):
    # Neither zero gain nor the finest available value is the continuum limit.
    h = np.array([0.01, 0.005, 0.0025, 0.00125])
    gain = limit + amplitude * (h / h[0]) ** order
    permutation = [2, 0, 3, 1]
    fitted = runner.fit_thermal_gain(
        h[permutation], gain[permutation], energy_scale=1e-3
    )
    assert fitted["status"] == "fit"
    assert fitted["order"] == pytest.approx(order, rel=1e-6)
    assert fitted["limit_J"] == pytest.approx(limit, rel=1e-6)
    assert fitted["amplitude_J"] == pytest.approx(amplitude, rel=1e-6)
    assert fitted["reference_spacing_m"] == h[0]
    assert fitted["residual_rms_J"] < 1e-13
    assert fitted["R2"] == pytest.approx(1)


@pytest.mark.parametrize(
    "gain,status",
    [
        ([1, 2], "at least three grids required"),
        ([1, 1, 1], "roundoff dominated"),
        ([1, 1 + 1e-15, 1 + 2e-15], "roundoff dominated"),
        ([1, 2, 1], "oscillatory gains"),
        ([1, 2, 4, 8], "no finite positive-order power-law fit"),
    ],
)
def test_thermal_gain_fit_reports_unavailable_sequences(runner, gain, status):
    h = 0.01 * 0.5 ** np.arange(len(gain))
    result = runner.fit_thermal_gain(h, gain)
    assert result["status"] == status
    assert "order" not in result
    assert "limit_J" not in result
    # Cancellation against a much larger initial inventory also limits accuracy.
    result = runner.fit_thermal_gain(
        [1, 0.5, 0.25], [1e-10, 2e-10, 3e-10], energy_scale=1e6
    )
    assert result["status"] == "roundoff dominated"


@pytest.mark.parametrize(
    "spacing,gain",
    [
        ([1, 0, 0.25], [1, 2, 3]),
        ([1, 1, 0.5], [1, 2, 3]),
        ([1, 0.5, 0.25], [1, np.nan, 3]),
    ],
)
def test_thermal_gain_fit_rejects_invalid_data(runner, spacing, gain):
    with pytest.raises(AuditError, match="distinct positive finite spacings"):
        runner.fit_thermal_gain(spacing, gain)


@pytest.mark.parametrize(
    "grids,expected",
    [
        ([8], None),
        ([8, 16], None),
        ([8, 16, 32], 2),
        ([8, 32, 128], 2),
        ([8, 16, 64], None),
    ],
)
def test_completed_subset_orders_use_actual_refinement(runner, grids, expected):
    rows = [(("fixed", "total_energy", nx, nx), {"Q": 1 + nx**-2}) for nx in grids]
    sequence = runner.grid_differences(rows, 2)
    if expected is None:
        assert sequence[0]["observed_order"] is None
    else:
        assert sequence[0]["observed_order"] == pytest.approx(expected)
    assert sequence[-1]["finest_difference_J"] is None
    if grids == [8, 16, 64]:
        assert sequence[0]["order_status"] == "unequal refinement ratios"


def test_relaxing_input_uses_formulation_boundary_parameter(runner):
    config = runner.configuration(runner.parser().parse_args([]))
    for formulation, expected, unused in (
        ("total_energy", "phi_boundary_velocity", "phi_boundary_timescale"),
        ("pressure", "phi_boundary_timescale", "phi_boundary_velocity"),
    ):
        text = runner.input_text(
            config,
            {"nx": 16, "nz": 16, "formulation": formulation, "boundary": "relax"},
        )
        assert expected + " = " in text
        assert unused + " = " not in text


@pytest.mark.parametrize(
    "boundary, radial_boundary", [("fixed", "fixed"), ("relax", "no_slip")]
)
def test_total_energy_input_uses_plan_c_boundary_options(
    runner, boundary, radial_boundary
):
    config = runner.configuration(
        runner.parser().parse_args(
            [
                "--levels",
                "2",
                "--base-nx",
                "16",
                "--base-nz",
                "16",
                "--phi-boundary",
                boundary,
                "--end-time",
                "10000",
                "--outputs",
                "50",
            ]
        )
    )
    for row in runner.grid_rows(config):
        text = runner.input_text(config, row)
        if row["formulation"] == "pressure":
            assert "radial_boundary =" not in text
            assert "inner_boundary_flags = 16" in text
            assert "outer_boundary_flags = 16" in text
            continue
        assert f"radial_boundary = {radial_boundary}\n" in text
        for removed in (
            "inner_boundary_flags",
            "outer_boundary_flags",
            "phi_boundary_relax",
            "phi_boundary_timescale",
            "phi_core_averagey",
            "neumann_boundary_average_z",
            "bndry_xin",
            "bndry_xout",
        ):
            assert removed not in text
        assert ("phi_boundary_velocity =" in text) == (radial_boundary == "no_slip")


@pytest.mark.parametrize("nx", [8, 16, 32, 64, 128, 256])
def test_tapered_initial_state_has_no_wall_sheet_or_total_energy_flow(runner, nx):
    from types import SimpleNamespace

    _, geometry = grid(nx=nx, nz=nx)
    run = SimpleNamespace(
        data=geometry,
        mxg=2,
        formulation="total_energy",
        scalar_option=lambda name: 1,
    )
    initial = runner.initial_state(run)
    interior = lambda field: field.isel(x=slice(2, -2))
    np.testing.assert_allclose(initial["phi"] + initial["Ph+"], 1, atol=1e-15)
    np.testing.assert_allclose(interior(initial["Vort"]), 0, atol=1e-12)
    np.testing.assert_allclose(
        interior(initial["Eion"] - 1.5 * initial["Ph+"]), 0, atol=1e-15
    )
    perturbation = interior(initial["Ph+"]) - 1
    # Both wall-adjacent rows lie in a constant buffer.
    assert float(abs(perturbation.isel(x=[0, -1])).max()) == 0
    run.formulation = "pressure"
    pressure_initial = runner.initial_state(run)
    if nx >= 16:
        # At these resolutions their stencil neighbours are also in the buffer.
        np.testing.assert_allclose(
            interior(pressure_initial["Vort"]).isel(x=[0, -1]), 0, atol=1e-12
        )
    # The legacy electron mass contribution leaves smooth bulk vorticity.
    expected = legacy_operator(1, initial["Pe"], geometry).divergence / -1836
    np.testing.assert_allclose(
        interior(pressure_initial["Vort"]), interior(expected), atol=1e-10
    )
    for name in ("Pe", "Ph+", "phi", "Nh+"):
        np.testing.assert_array_equal(
            interior(initial[name]), interior(pressure_initial[name])
        )


def test_initial_pressure_preserves_central_perturbation_and_flat_walls(runner):
    x = np.array([0, 0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 0.8, 0.9, 0.95, 1])
    z = np.arange(16) * 2 * np.pi / 16
    pressure = runner.initial_pressure(x, z)
    np.testing.assert_array_equal(pressure[[0, 1, 2, 8, 9, 10]], 1)
    original = 1 + 0.1 * np.sin(2 * np.pi * x[:, None]) * np.sin(
        2 * np.pi * x[:, None] - z[None, :]
    )
    np.testing.assert_allclose(pressure[3:8], original[3:8], atol=1e-15)


def test_saved_time_lookup_rejects_nearest_output(runner):
    with pytest.raises(AuditError):
        runner.exact_index(np.arange(5), 2.2)
    assert runner.exact_index(np.arange(5), 2) == 2


def test_singleton_y_fluxes_have_finite_zero_divergence():
    from energy_analysis.transport import flux_dataset, parallel_flux

    field, geometry = grid()
    geometry["dy"] = xr.ones_like(field)
    geometry["g_22"] = xr.ones_like(field)
    result = parallel_flux(field + 1, xr.zeros_like(field), field + 1, geometry, myg=0)
    assert np.isfinite(result.divergence).all()
    assert np.max(abs(result.divergence)) == 0
    result = flux_dataset({"y": xr.zeros_like(field)}, geometry)
    assert np.max(abs(result.divergence)) == 0


def test_run_timeout_records_failure_and_stops_process(runner, tmp_path):
    import json

    exe = tmp_path / "delayed-simulation"
    exe.write_text("#!/usr/bin/env python3\nimport time\ntime.sleep(30)\n")
    exe.chmod(0o755)
    directory = tmp_path / "run with spaces"
    directory.mkdir()
    with pytest.raises(AuditError, match="Run timed out after 0.1 s"):
        runner.execute(exe, directory, {"nproc": 1}, 0.1)
    status = json.loads((directory / "status.json").read_text())
    assert status["exit_code"] == -1
    assert status["termination_reason"] == "timeout"
    assert status["command"] == [str(exe), "-d", str(directory)]
    assert status["elapsed_s"] < 6


def test_native_library_changes_are_part_of_runtime_hash(runner, tmp_path):
    exe = tmp_path / "hermes-3"
    exe.write_bytes(b"simulation")
    library = tmp_path / "libhermes-3-lib.so"
    library.write_bytes(b"first")
    before = runner.runtime_artifacts(exe)
    library.write_bytes(b"second")
    after = runner.runtime_artifacts(exe)
    assert before[str(exe)] == after[str(exe)]
    assert before[str(library)] != after[str(library)]


def test_new_study_uses_shared_analysis_without_xhermes(runner, tmp_path, monkeypatch):
    import json

    monkeypatch.setitem(sys.modules, "xhermes", None)
    executable = tmp_path / "hermes-3"
    executable.write_bytes(b"simulation")
    executable.chmod(0o755)

    def startup(exe, directory, config, timeout):
        (directory / "initial-state.nc").write_bytes(b"matched state")
        for rank in range(config["nproc"]):
            (directory / f"BOUT.restart.{rank}.nc").write_bytes(b"restart")

    monkeypatch.setattr(runner, "ensure_startup", startup)
    output = tmp_path / "results"
    assert (
        runner.main(
            [
                "--exe",
                str(executable),
                "--output-root",
                str(output),
                "--smoke",
                "--prepare-only",
                "--nproc",
                "1",
            ]
        )
        == 0
    )
    study = next(output.glob("interchange-*"))
    manifest = json.loads((study / "manifest.json").read_text())
    assert "energy_analysis" in manifest["revisions"]
    assert "xhermes" not in manifest["environment"]
    for case in manifest["runs"]:
        assert "${" not in (study / case["directory"] / "BOUT.inp").read_text()
        assert (study / case["directory"] / "remote-case.json").is_file()


def test_reanalysis_preserves_saved_style_and_allows_overrides(runner, tmp_path):
    import json

    import matplotlib

    config = {
        "field": "ion_pressure",
        "snapshot_times": [0, 3, 6, 12],
        "convergence_time": 6,
    }
    style = tmp_path / "style.json"
    style.write_text(json.dumps({"font.size": 13, "axes.labelsize": 11}))
    saved = tmp_path / "figure-config.json"
    with matplotlib.rc_context():
        runner.figures(tmp_path, config, {}, style)
        original = json.loads(saved.read_text())
        runner.figures(tmp_path, config, {})
        assert json.loads(saved.read_text()) == original
        style.write_text(json.dumps({"font.size": 10}))
        runner.figures(tmp_path, config, {}, style)
        updated = json.loads(saved.read_text())["rcParams"]
        assert updated["font.size"] == 10
        assert updated["axes.labelsize"] == 11


@pytest.mark.parametrize("energy_format", ["png", "pdf-panels"])
def test_energy_format_is_plotting_only_and_reaches_reanalysis(
    runner, tmp_path, monkeypatch, energy_format
):
    import json

    default = runner.parser().parse_args([])
    assert default.energy_format == "png"
    argv = ["--analyse-only", str(tmp_path), "--energy-format", energy_format]
    args = runner.parser().parse_args(argv)
    assert runner.configuration(args) == runner.configuration(default)
    (tmp_path / "manifest.json").write_text(
        json.dumps({"schema_version": runner.SCHEMA})
    )
    calls = []
    monkeypatch.setattr(runner, "analyse", lambda *a, **kw: calls.append(kw))
    assert runner.main(argv) == 0
    assert calls[0]["energy_format"] == energy_format
    assert calls[0]["allow_incomplete"] is True


def synthetic_plot_result(nx, gain):
    """Inventories with a larger intermediate change than the final signed gain."""
    inventory = xr.Dataset(
        {
            "ion_thermal": ("t", [10, 20, 15, 10 + gain]),
            "electron_thermal": ("t", [5, 5, 5, 5]),
            "ion_perpendicular": ("t", [0, 1, 2, 3]),
            "plasma_reconstructed": ("t", [15, 26, 22, 18 + gain]),
            "plasma_native": ("t", [15, 25, 21, 17 + gain]),
        },
        coords={"t": [0, 0.001, 0.002, 0.004]},
    )
    return {
        "inventory": inventory,
        "lx": 0.05,
        "lz": 0.05,
        "x": (np.arange(nx) + 0.5) * 0.05 / nx,
        "z": np.arange(nx) * 0.05 / nx,
        "snapshots": [np.ones((nx, nx)) * p for p in (1, 2, 3, 4)],
        "kinetic_snapshots": [np.ones((nx, nx)) * p for p in (0, 1, 2, 3)],
        "comparison": np.ones((nx, nx)),
        "volume": np.ones((nx, nx)),
        "unit": "Pa",
        "Q": 1 + nx**-2,
        "report": {
            "volume_m3": 0.001,
            "ion_thermal_gain_J": gain,
            "final_drift_J": 2 + gain,
            "max_relative_drift": 0.5,
            "final_corrected_residual_J": None,
            "corrected_residual_J": None,
        },
    }


def test_thermal_gain_pdf_overlays_signed_gains_and_continuum_fit(
    runner, tmp_path, monkeypatch
):
    import json

    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure

    groups = {
        "pressure": [
            (("fixed", "pressure", nx, nx), synthetic_plot_result(nx, -8 / nx))
            for nx in (8, 16)
        ],
        "total_energy": [
            (
                ("fixed", "total_energy", nx, nx),
                synthetic_plot_result(nx, 0.25 + (8 / nx) ** 2),
            )
            for nx in (8, 16, 32)
        ],
    }
    captured = []
    original_savefig = Figure.savefig

    def capture(fig, path, **kw):
        captured.append(fig.axes[0])
        original_savefig(fig, path, **kw)

    monkeypatch.setattr(Figure, "savefig", capture)
    with plt.rc_context({"figure.constrained_layout.use": True, "pdf.fonttype": 42}):
        runner.thermal_gain_figure(tmp_path, {"end_time": 4}, "fixed", groups)
    assert (tmp_path / "ion-thermal-gain-fixed.pdf").read_bytes().startswith(b"%PDF-")
    (axis,) = captured
    assert axis.get_xscale() == "log"
    assert axis.get_yscale() == "linear"
    lines = {line.get_label(): line for line in axis.lines}
    np.testing.assert_allclose(lines["Pressure evolution"].get_xdata(), [6.25, 3.125])
    np.testing.assert_allclose(lines["Pressure evolution"].get_ydata(), [-1, -0.5])
    np.testing.assert_allclose(
        lines["Total ion energy"].get_ydata(), [1.25, 0.5, 0.3125]
    )
    fit_line = next(
        line for label, line in lines.items() if label.startswith("Ion-energy fit")
    )
    assert fit_line.get_xdata().min() < 0.05 / 32 * 1e3
    limit_line = next(line for label, line in lines.items() if "h\\to0" in label)
    np.testing.assert_allclose(limit_line.get_ydata(), 0.25, atol=1e-7)
    report = json.loads((tmp_path / "ion-thermal-gain-fixed.json").read_text())
    assert report["ion_energy_fit"]["order"] == pytest.approx(2, rel=1e-6)
    assert report["ion_energy_fit"]["limit_J"] == pytest.approx(0.25, rel=1e-6)
    assert report["data"][0]["ion_thermal_gain_J"] == -1
    assert report["data"][0]["initial_ion_thermal_J"] == 10
    assert report["data"][0]["final_ion_thermal_J"] == 9


@pytest.mark.parametrize("formulation", ["pressure", "total_energy"])
def test_partial_reanalysis_plots_single_completed_formulation(
    runner, tmp_path, monkeypatch, formulation
):
    import json

    config = runner.configuration(runner.parser().parse_args(["--smoke"]))
    rows = list(runner.grid_rows(config))
    complete = next(row for row in rows if row["formulation"] == formulation)
    for row in rows:
        directory = tmp_path / row["directory"]
        directory.mkdir()
        (directory / "BOUT.inp").write_text(runner.input_text(config, row))
        row["input_sha256"] = runner.digest(directory / "BOUT.inp")
    runner.save_json(tmp_path / complete["directory"] / "status.json", {"exit_code": 0})
    monkeypatch.setattr(
        runner, "analyse_run", lambda *args: synthetic_plot_result(8, -0.5)
    )
    manifest = {"configuration": config, "runs": rows}
    results = runner.analyse(tmp_path, manifest, allow_incomplete=True)
    assert len(results) == 1
    boundary = complete["boundary"]
    report = json.loads((tmp_path / f"ion-thermal-gain-{boundary}.json").read_text())
    assert [sample["formulation"] for sample in report["data"]] == [formulation]
    assert report["data"][0]["ion_thermal_gain_J"] == -0.5
    assert report["ion_energy_fit"]["status"] == (
        "at least three grids required"
        if formulation == "total_energy"
        else "no total-ion-energy runs available"
    )
    convergence = json.loads((tmp_path / "convergence-table.json").read_text())
    assert convergence[0]["observed_order"] is None
    assert convergence[0]["L2"] is None
    assert (tmp_path / f"state-{boundary}.pdf").is_file()
    assert (tmp_path / f"ion-thermal-gain-{boundary}.pdf").is_file()
    # Partial analysis must retain the input hash guard, including pending runs.
    pending = next(row for row in rows if row is not complete)
    (tmp_path / pending["directory"] / "BOUT.inp").write_text("changed")
    with pytest.raises(AuditError, match="Input was modified"):
        runner.analyse(tmp_path, manifest, allow_incomplete=True)


def test_partial_reanalysis_requires_a_completed_run(runner, tmp_path):
    config = runner.configuration(runner.parser().parse_args(["--smoke"]))
    row = next(runner.grid_rows(config))
    directory = tmp_path / row["directory"]
    directory.mkdir()
    (directory / "BOUT.inp").write_text(runner.input_text(config, row))
    row["input_sha256"] = runner.digest(directory / "BOUT.inp")
    with pytest.raises(AuditError, match="No completed runs"):
        runner.analyse(
            tmp_path, {"configuration": config, "runs": [row]}, allow_incomplete=True
        )


def test_energy_exports_keep_panel_data_and_standalone_legends(
    runner, tmp_path, monkeypatch
):
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure

    groups = {}
    for form in ("pressure", "total_energy"):
        groups[form] = []
        for nx in (8, 16):
            inv = xr.Dataset(
                {
                    "ion_thermal": ("t", [10, 11, 13]),
                    "electron_thermal": ("t", [5, 4, 2]),
                    "ion_perpendicular": ("t", [0, 2, 4]),
                    "plasma_reconstructed": ("t", [15, 17, 19]),
                    "plasma_native": ("t", [15, 16, 18]),
                },
                coords={"t": [0, 0.001, 0.002]},
            )
            groups[form].append(
                (
                    ("fixed", form, nx, nx),
                    {
                        "inventory": inv,
                        "report": {
                            "corrected_residual_J": [0, 0.3, 0.6]
                            if form == "total_energy"
                            else None
                        },
                    },
                )
            )

    captured = {}
    original_savefig = Figure.savefig

    def capture(fig, path, **kw):
        original_savefig(fig, path, **kw)
        captured[Path(path).name] = [
            {
                "title": axis.get_title(),
                "xlabel": axis.get_xlabel(),
                "legend": axis.get_legend(),
                "lines": [
                    (
                        np.asarray(line.get_xdata()).copy(),
                        np.asarray(line.get_ydata()).copy(),
                        line.get_color(),
                        line.get_linestyle(),
                        line.get_label(),
                    )
                    for line in axis.lines
                ],
            }
            for axis in fig.axes
        ]

    monkeypatch.setattr(Figure, "savefig", capture)
    combined, separate = tmp_path / "png", tmp_path / "pdf"
    combined.mkdir()
    separate.mkdir()
    with plt.rc_context({"figure.constrained_layout.use": True}):
        runner.energy_figures(combined, {"viscosity": 1}, "fixed", groups)
        runner.energy_figures(
            separate, {"viscosity": 1}, "fixed", groups, energy_format="pdf-panels"
        )

    assert [p.name for p in combined.iterdir()] == ["energy-fixed.png"]
    assert (combined / "energy-fixed.png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    names = (
        "ion-thermal",
        "electron-thermal",
        "perpendicular-kinetic",
        "shared-relative-change",
        "model-relative-change",
        "estimator-difference",
    )
    assert {p.name for p in separate.iterdir()} == {
        f"energy-fixed-{name}.pdf" for name in names
    }
    expected = (
        [0, 1, 3],
        [0, -1, -3],
        [0, 2, 4],
        [0, 2 / 15, 4 / 15],
        [0, 1 / 15, 3 / 15],
        [0, 1, 1],
    )
    for panel, name, values in zip(captured["energy-fixed.png"], names, expected):
        path = separate / f"energy-fixed-{name}.pdf"
        assert path.read_bytes().startswith(b"%PDF-")
        (standalone,) = captured[path.name]
        assert standalone["title"] == panel["title"]
        assert standalone["xlabel"] == panel["xlabel"] == "Time [ms]"
        assert standalone["legend"] is not None
        assert (
            len(standalone["lines"])
            == len(panel["lines"])
            == (6 if name == "model-relative-change" else 4)
        )
        np.testing.assert_allclose(standalone["lines"][0][0], [0, 1, 2])
        np.testing.assert_allclose(standalone["lines"][0][1], values)
        for split_line, combined_line in zip(standalone["lines"], panel["lines"]):
            np.testing.assert_array_equal(split_line[0], combined_line[0])
            np.testing.assert_array_equal(split_line[1], combined_line[1])
            assert split_line[2:] == combined_line[2:]
        if name == "model-relative-change":
            np.testing.assert_allclose(standalone["lines"][-1][1], [0, 0.02, 0.04])


def write_publication_dump(directory, rank, times, *, bad_final=False, seed=False):
    fields = ("Ne", "Nh+", "Pe", "phi", "Vort", "Ph+")
    values = np.ones((len(times), 8, 1, 4))
    # Physical state is required; undefined diagnostic/guard entries are allowed.
    values[:, :2] = values[:, -2:] = np.nan
    ds = xr.Dataset(
        {name: (("t", "x", "y", "z"), values.copy()) for name in fields},
        coords={"t": times},
    )
    if bad_final:
        ds["Ph+"][-1, 3, 0, 1] = np.nan
    ds.to_netcdf(directory / f"BOUT.dmp.{rank}.nc")
    if seed:
        xr.Dataset({"tt": 0.0}).to_netcdf(directory / f"BOUT.restart.{rank}.nc")


@pytest.fixture
def completed_output_case(tmp_path):
    directory = tmp_path / "run"
    directory.mkdir()
    config = {"nproc": 2, "end_time": 4, "outputs": 4}
    for rank in range(2):
        write_publication_dump(directory, rank, np.arange(5.0))
    return directory, config


def test_nonzero_exit_accepts_complete_outputs(
    runner, completed_output_case, monkeypatch
):
    import json

    directory, config = completed_output_case
    # Use a real process return code without requiring MPI for this test.
    executable = directory / "finished-simulation"
    executable.write_text("#!/bin/sh\nexit 143\n")
    executable.chmod(0o755)
    monkeypatch.setitem(config, "nproc", 1)
    (directory / "BOUT.dmp.1.nc").unlink()
    status = runner.execute(executable, directory, config, 2)
    assert status["exit_code"] == 143
    assert status["outputs_complete"]
    assert status["output_check"] == {
        "ranks": 1,
        "records_per_rank": 5,
        "final_time_normalised": 4,
    }
    assert json.loads((directory / "status.json").read_text()) == status


@pytest.mark.parametrize("defect", ["rank", "time", "field", "unreadable"])
def test_complete_output_check_rejects_incomplete_ranks(
    runner, completed_output_case, defect
):
    import json

    directory, config = completed_output_case
    file = directory / "BOUT.dmp.1.nc"
    if defect == "rank":
        file.unlink()
    elif defect == "time":
        write_publication_dump(directory, 1, np.arange(4.0))
    elif defect == "field":
        write_publication_dump(directory, 1, np.arange(5.0), bad_final=True)
    else:
        file.write_bytes(b"unreadable output")
    runner.save_json(directory / "status.json", {"exit_code": 143})
    assert not runner.completed_run(directory, config)
    status = json.loads((directory / "status.json").read_text())
    assert status["exit_code"] == 143
    assert not status["outputs_complete"]
    assert status["output_check_error"]


def test_resume_recovers_old_status_and_rechecks_saved_completion(
    runner, completed_output_case
):
    directory, config = completed_output_case
    runner.save_json(directory / "status.json", {"exit_code": 143})
    assert runner.completed_run(directory, config)
    write_publication_dump(directory, 1, np.arange(4.0))
    assert not runner.completed_run(directory, config)


def test_seed_completion_requires_restart_files(runner, completed_output_case):
    directory, config = completed_output_case
    for rank in range(2):
        write_publication_dump(directory, rank, [0.0], seed=True)
    assert (
        runner.verify_completed_outputs(directory, config, seed=True)[
            "records_per_rank"
        ]
        == 1
    )
    (directory / "BOUT.restart.1.nc").unlink()
    with pytest.raises(AuditError, match="restart per MPI rank"):
        runner.verify_completed_outputs(directory, config, seed=True)


def test_interrupt_does_not_accept_completed_outputs(runner, tmp_path, monkeypatch):
    class Process:
        pid = 123
        calls = 0

        def wait(self, timeout):
            self.calls += 1
            if self.calls == 1:
                raise KeyboardInterrupt
            return 0

    monkeypatch.setattr(runner.subprocess, "Popen", lambda *a, **k: Process())
    monkeypatch.setattr(runner.os, "killpg", lambda *a: None)

    def unexpected_check(*args, **kwargs):
        pytest.fail("Interrupt must stop the invocation before completion recovery")

    monkeypatch.setattr(runner, "accept_completed_outputs", unexpected_check)
    with pytest.raises(KeyboardInterrupt):
        runner.execute(tmp_path / "simulation", tmp_path, {"nproc": 1}, 1)


@pytest.fixture
def resume_case(runner, tmp_path, monkeypatch):
    executable = tmp_path / "hermes-3"
    executable.write_bytes(b"simulation")
    study = tmp_path / "study"
    study.mkdir()
    argv = ["--exe", str(executable), "--resume", str(study)]
    manifest = {
        "configuration": runner.configuration(runner.parser().parse_args(argv)),
        "executable": {"sha256": runner.digest(executable)},
        "runtime_artifacts": runner.runtime_artifacts(executable),
        "template_sha256": runner.digest(
            runner.HERE / "publication-input.inp.template"
        ),
        "runner_sha256": "original-runner",
        "runs": [],
    }
    runner.save_json(study / "manifest.json", manifest)
    monkeypatch.setattr(runner, "analyse", lambda *a, **k: {})
    return study, argv, manifest


def test_resume_runner_change_is_explicit_and_recorded(runner, resume_case):
    import json

    study, argv, manifest = resume_case
    with pytest.raises(AuditError, match="Runner changed"):
        runner.main(argv)
    assert runner.main([*argv, "--allow-runner-change"]) == 0
    updated = json.loads((study / "manifest.json").read_text())
    assert updated["runner_sha256"] == manifest["runner_sha256"]
    assert updated["runner_updates"][0]["previous_sha256"] == "original-runner"
    assert updated["resume_runner_sha256"] == runner.digest(runner.__file__)
    assert runner.main(argv) == 0
    assert json.loads((study / "manifest.json").read_text()) == updated


@pytest.mark.parametrize(
    "change", ["configuration", "executable", "template", "library"]
)
def test_runner_override_preserves_simulation_hash_guards(runner, resume_case, change):
    study, argv, manifest = resume_case
    if change == "configuration":
        manifest["configuration"]["end_time"] *= 2
    elif change == "executable":
        manifest["executable"]["sha256"] = "changed"
    elif change == "template":
        manifest["template_sha256"] = "changed"
    else:
        manifest["runtime_artifacts"] = {}
    runner.save_json(study / "manifest.json", manifest)
    with pytest.raises(
        AuditError, match="matching configuration, executable and template"
    ):
        runner.main([*argv, "--allow-runner-change"])


@pytest.fixture
def portable_study(runner, tmp_path):
    exe = tmp_path / "hermes-3"
    exe.write_bytes(b"simulation")
    study = tmp_path / "study"
    study.mkdir()
    config = runner.configuration(
        runner.parser().parse_args(
            [
                "--levels",
                "2",
                "--base-nx",
                "4",
                "--base-nz",
                "4",
                "--end-time",
                "4",
                "--outputs",
                "4",
            ]
        )
    )
    rows = list(runner.grid_rows(config))
    runner.add_inputs(study, config, rows)
    manifest = {
        "schema_version": runner.SCHEMA,
        "created_utc": "test-study",
        "configuration": config,
        "executable": {"path": str(exe), "sha256": runner.digest(exe)},
        "runtime_artifacts": runner.runtime_artifacts(exe),
        "template_sha256": runner.digest(
            runner.HERE / "publication-input.inp.template"
        ),
        "runner_sha256": runner.digest(runner.__file__),
        "revisions": {"hermes": {"revision": "expected-source"}},
        "runs": rows,
    }
    runner.save_json(study / "manifest.json", manifest)
    return study, exe, manifest


def test_extension_runs_only_new_grids_and_keeps_completed_data(
    runner, portable_study, monkeypatch
):
    import json

    study, exe, manifest = portable_study
    original = {}
    for row in manifest["runs"]:
        directory = study / row["directory"]
        runner.save_json(directory / "status.json", {"exit_code": 0})
        (directory / "inventory.nc").write_bytes(b"original analysis")
        for name in ("BOUT.inp", "status.json", "inventory.nc"):
            original[directory / name] = (directory / name).read_bytes()
    prepared, executed, analyses = [], [], []
    monkeypatch.setattr(
        runner, "ensure_startup", lambda e, d, c, t: prepared.append(d.name)
    )

    def execute(e, d, c, t):
        executed.append(d.name)
        runner.save_json(d / "status.json", {"exit_code": 0})

    monkeypatch.setattr(runner, "execute", execute)
    monkeypatch.setattr(runner, "analyse", lambda s, m, **k: analyses.append(m["runs"]))
    assert runner.main(["--extend", str(study), "--levels", "3"]) == 0
    assert set(executed) == {"fixed-pressure-16x16", "fixed-total_energy-16x16"}
    assert prepared == executed
    updated = json.loads((study / "manifest.json").read_text())
    assert updated["configuration"]["levels"] == 3
    assert len(updated["runs"]) == len(analyses[0]) == 6
    assert updated["extensions"][0]["previous_levels"] == 2
    for path, content in original.items():
        assert path.read_bytes() == content
    # Ordinary resume also accepts a larger level count with matching flags.
    flags = [
        "--resume",
        str(study),
        "--exe",
        str(exe),
        "--levels",
        "4",
        "--base-nx",
        "4",
        "--base-nz",
        "4",
        "--end-time",
        "4",
        "--outputs",
        "4",
    ]
    executed.clear()
    assert runner.main(flags) == 0
    assert set(executed) == {"fixed-pressure-32x32", "fixed-total_energy-32x32"}


@pytest.mark.parametrize(
    "flags", [["--levels", "1"], ["--levels", "3", "--viscosity", "2"]]
)
def test_extension_rejects_shrinking_and_other_setting_changes(
    runner, portable_study, flags
):
    study, _, _ = portable_study
    original = (study / "manifest.json").read_bytes()
    with pytest.raises(AuditError):
        runner.main(["--extend", str(study), *flags])
    assert (study / "manifest.json").read_bytes() == original


def test_extension_does_not_adopt_conflicting_input(runner, portable_study):
    study, _, manifest = portable_study
    directory = study / "fixed-total_energy-16x16"
    directory.mkdir()
    (directory / "BOUT.inp").write_text("different experiment")
    original = (study / "manifest.json").read_bytes()
    with pytest.raises(AuditError, match="Existing input differs"):
        runner.extend_study(study, manifest, {**manifest["configuration"], "levels": 3})
    assert (study / "manifest.json").read_bytes() == original
    assert not (study / "fixed-pressure-16x16").exists()


def test_prepare_selection_preserves_restarts_for_resume(
    runner, portable_study, monkeypatch
):
    study, exe, manifest = portable_study
    calls = []

    def startup(e, d, c, t):
        calls.append(d.name)
        (d / "startup").mkdir()
        xr.Dataset({"tt": 0.0}).to_netcdf(d / "BOUT.restart.0.nc")
        xr.Dataset({"phi": 0.0}).to_netcdf(d / "initial-state.nc")

    monkeypatch.setattr(runner, "build_startup", startup)
    monkeypatch.setattr(
        runner,
        "execute",
        lambda *a, **k: pytest.fail("Full simulation ran during preparation"),
    )
    monkeypatch.setattr(
        runner, "analyse", lambda *a, **k: pytest.fail("Incomplete study was analysed")
    )
    name = "fixed-pressure-16x16"
    flags = ["--extend", str(study), "--levels", "3", "--prepare-only", "--cases", name]
    assert runner.main(flags) == 0
    assert calls == [name]
    directory = study / name
    original = (directory / "BOUT.restart.0.nc").read_bytes()
    runner.ensure_startup(exe, directory, {**manifest["configuration"], "levels": 3}, 2)
    assert calls == [name]
    assert (directory / "BOUT.restart.0.nc").read_bytes() == original
    assert (directory / "run-remote.py").is_file()
    assert not list(directory.glob("failed-*"))
    (directory / "BOUT.dmp.0.nc").write_bytes(b"returned output")
    with pytest.raises(AuditError, match="Unregistered outputs"):
        runner.ensure_startup(exe, directory, manifest["configuration"], 2)


@pytest.fixture
def importable_case(runner, portable_study, monkeypatch, tmp_path):
    from dataclasses import dataclass

    study, _, manifest = portable_study
    row = manifest["runs"][0]
    destination = study / row["directory"]
    returned = tmp_path / "returned"
    source = returned / row["directory"]
    source.mkdir(parents=True)
    (source / "BOUT.inp").write_bytes((destination / "BOUT.inp").read_bytes())
    write_publication_dump(source, 0, np.arange(5.0))
    with xr.open_dataset(source / "BOUT.dmp.0.nc") as ds:
        initial = ds.isel(t=0, drop=True).load()
    initial.to_netcdf(destination / "initial-state.nc")

    @dataclass
    class Run:
        data: xr.Dataset
        source: str
        formulation: str = "pressure"
        length_scale: float = 1

        def interior(self, value):
            return value.isel(x=slice(2, -2))

    def open_run(directory):
        from boutdata.data import BoutOptionsFile

        ds = xr.open_dataset(directory / "BOUT.dmp.0.nc")
        options = directory / "BOUT.settings"
        if not options.exists():
            options = directory / "BOUT.inp"
        ds.attrs["options"] = BoutOptionsFile(str(options))
        for name, value in {
            "dx": 0.05 / 4,
            "dz": 0.05 / 4,
            "g11": 1.0,
            "g33": 1.0,
        }.items():
            ds[name] = value
        return Run(ds, str(directory))

    monkeypatch.setattr(runner, "open_run", open_run)
    return study, manifest, row, returned, source


def test_import_bare_outputs_without_executable_and_resume_reuses_case(
    runner, importable_case, monkeypatch, capsys
):
    import json

    study, manifest, row, returned, source = importable_case
    (source / "BOUT.log.0").write_text("remote run")
    assert (
        runner.main(
            [
                "--import-runs",
                str(study),
                "--import-from",
                str(returned),
                "--cases",
                row["directory"],
                "--exe",
                "/missing/executable",
            ]
        )
        == 0
    )
    assert "provenance are unknown" in capsys.readouterr().err
    destination = study / row["directory"]
    status = json.loads((destination / "status.json").read_text())
    assert status["exit_code"] is None
    assert status["origin"] == "external"
    assert status["outputs_complete"]
    assert status["initial_state_check"]["ion_pressure"] == 0
    assert (destination / "BOUT.log.0").read_text() == "remote run"
    assert runner.completed_run(destination, manifest["configuration"])
    assert runner.import_runs(study, manifest, [row["directory"]], returned) == 0
    # Imported files are validated again on resume, even with unknown exit code.
    write_publication_dump(destination, 0, np.arange(4.0))
    with pytest.raises(AuditError, match="output times"):
        runner.completed_run(destination, manifest["configuration"])


@pytest.mark.parametrize(
    "defect", ["time", "initial", "nan-initial", "input", "rank", "grid"]
)
def test_import_rejects_invalid_outputs_before_registering(
    runner, importable_case, defect
):
    study, manifest, row, returned, source = importable_case
    file = source / "BOUT.dmp.0.nc"
    if defect == "time":
        write_publication_dump(source, 0, np.arange(4.0))
    elif defect == "input":
        (source / "BOUT.inp").write_text("different settings")
    elif defect == "rank":
        write_publication_dump(source, 1, np.arange(5.0))
    else:
        with xr.open_dataset(file) as ds:
            ds = ds.load()
        if defect == "grid":
            ds = ds.isel(z=slice(0, 3))
        else:
            ds["phi"][0, 3, 0, 1] = np.nan if defect == "nan-initial" else 2
        ds.to_netcdf(file)
    with pytest.raises(AuditError):
        runner.import_runs(study, manifest, [row["directory"]], returned)
    assert not (study / row["directory"] / "status.json").exists()
    assert not list(study.glob(".import-*"))
    assert not (study / row["directory"] / "BOUT.dmp.0.nc").exists()


def test_import_files_copied_directly_back(runner, importable_case):
    import shutil

    study, manifest, row, _, source = importable_case
    destination = study / row["directory"]
    shutil.copy2(source / "BOUT.dmp.0.nc", destination / "BOUT.dmp.0.nc")
    assert runner.import_runs(study, manifest, [row["directory"]]) == 1
    assert runner.completed_run(destination, manifest["configuration"])


def test_import_archives_failed_outputs_logs_and_stale_settings(
    runner, importable_case
):
    import json

    study, manifest, row, returned, _ = importable_case
    destination = study / row["directory"]
    write_publication_dump(destination, 0, np.arange(4.0))
    runner.save_json(destination / "status.json", {"exit_code": -1})
    (destination / "BOUT.settings").write_text("old settings")
    (destination / "run.log").write_text("failed attempt")
    assert runner.import_runs(study, manifest, [row["directory"]], returned) == 1
    (archive,) = destination.glob("pre-import-*")
    assert (archive / "BOUT.dmp.0.nc").is_file()
    assert (archive / "BOUT.settings").read_text() == "old settings"
    assert (archive / "run.log").read_text() == "failed attempt"
    assert not (destination / "BOUT.settings").exists()
    status = json.loads((destination / "status.json").read_text())
    assert status["previous_status"]["exit_code"] == -1
    assert status["options_source"] == "BOUT.inp"


def test_import_rejects_command_line_override_in_returned_settings(
    runner, importable_case
):
    study, manifest, row, returned, source = importable_case
    (source / "BOUT.settings").write_text(
        (source / "BOUT.inp").read_text().replace("viscosity = 1", "viscosity = 2")
    )
    with pytest.raises(AuditError, match="Returned setting differs"):
        runner.import_runs(study, manifest, [row["directory"]], returned)
    assert not (study / row["directory"] / "status.json").exists()


@pytest.fixture
def remote_runner(runner):
    spec = importlib.util.spec_from_file_location(
        "publication_remote", runner.HERE / "publication_remote.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_remote_runner_records_real_exit_and_build_provenance(
    runner, remote_runner, importable_case, tmp_path, monkeypatch
):
    import json
    import shutil

    study, manifest, row, returned, source = importable_case
    destination = study / row["directory"]
    xr.Dataset({"tt": 0.0}).to_netcdf(destination / "BOUT.restart.0.nc")
    runner.prepare_remote(destination, row, manifest)
    execution_dir = tmp_path / "remote-execution" / row["directory"]
    shutil.copytree(destination, execution_dir)
    exe = tmp_path / "finished-simulation"
    exe.write_text("#!/bin/sh\nexit 143\n")
    exe.chmod(0o755)
    monkeypatch.setattr(
        remote_runner,
        "revision",
        lambda p: {"revision": "expected-source", "dirty": False},
    )
    assert (
        remote_runner.main(
            ["--exe", str(exe), "--source-root", str(tmp_path)], directory=execution_dir
        )
        == 143
    )
    execution = json.loads((execution_dir / "remote-execution.json").read_text())
    assert execution["exit_code"] == 143
    assert (
        execution["input_files"]
        == json.loads((destination / "remote-case.json").read_text())["input_files"]
    )
    assert execution["runtime_artifacts"][str(exe)] == runner.digest(exe)
    shutil.copy2(
        execution_dir / "remote-execution.json", source / "remote-execution.json"
    )
    assert runner.import_runs(study, manifest, [row["directory"]], returned) == 1
    status = json.loads((destination / "status.json").read_text())
    assert status["exit_code"] == 143
    assert status["remote_execution"]["hostname"]
    with pytest.raises(ValueError, match="Outputs already exist"):
        remote_runner.main(["--exe", str(exe)], directory=execution_dir)


@pytest.mark.parametrize("defect", ["preparation", "revision", "runtime"])
def test_import_rejects_mismatched_execution_metadata(runner, importable_case, defect):
    import json

    study, manifest, row, returned, source = importable_case
    destination = study / row["directory"]
    xr.Dataset({"tt": 0.0}).to_netcdf(destination / "BOUT.restart.0.nc")
    runner.prepare_remote(destination, row, manifest)
    execution = {
        "schema_version": 1,
        "case": row["directory"],
        "exit_code": 0,
        "prepared_case_sha256": runner.digest(destination / "remote-case.json"),
        "input_files": json.loads((destination / "remote-case.json").read_text())[
            "input_files"
        ],
        "runtime_artifacts": {"exe": "build"},
        "runtime_artifacts_after": {"exe": "build"},
        "revisions": {"hermes": {"revision": "expected-source"}},
    }
    if defect == "preparation":
        execution["prepared_case_sha256"] = "other-case"
    elif defect == "revision":
        execution["revisions"]["hermes"]["revision"] = "different-source"
    else:
        execution["runtime_artifacts_after"] = {"exe": "rebuilt"}
    runner.save_json(source / "remote-execution.json", execution)
    with pytest.raises(AuditError):
        runner.import_runs(study, manifest, [row["directory"]], returned)
    assert not (destination / "status.json").exists()


def test_remote_helper_runs_without_site_packages(runner, portable_study, tmp_path):
    import json
    import subprocess
    import sys

    study, _, manifest = portable_study
    row = manifest["runs"][0]
    directory = study / row["directory"]
    xr.Dataset({"tt": 0.0}).to_netcdf(directory / "BOUT.restart.0.nc")
    xr.Dataset({"phi": 0.0}).to_netcdf(directory / "initial-state.nc")
    runner.prepare_remote(directory, row, manifest)
    exe = tmp_path / "remote-simulation"
    exe.write_text("#!/bin/sh\nexit 143\n")
    exe.chmod(0o755)
    result = subprocess.run(
        [
            sys.executable,
            "-I",
            "-S",
            str(directory / "run-remote.py"),
            "--exe",
            str(exe),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 143, result.stderr
    assert (
        json.loads((directory / "remote-execution.json").read_text())["exit_code"]
        == 143
    )


def test_remote_helper_rejects_known_source_change_before_launch(
    runner, remote_runner, portable_study, monkeypatch
):
    study, exe, manifest = portable_study
    row = manifest["runs"][0]
    directory = study / row["directory"]
    xr.Dataset({"tt": 0.0}).to_netcdf(directory / "BOUT.restart.0.nc")
    xr.Dataset({"phi": 0.0}).to_netcdf(directory / "initial-state.nc")
    runner.prepare_remote(directory, row, manifest)
    monkeypatch.setattr(
        remote_runner,
        "revision",
        lambda p: {"revision": "wrong-source", "dirty": False},
    )
    monkeypatch.setattr(
        remote_runner.subprocess,
        "Popen",
        lambda *a, **k: pytest.fail("Mismatched remote source was launched"),
    )
    with pytest.raises(ValueError, match="revision differs"):
        remote_runner.main(
            ["--exe", str(exe), "--source-root", str(study)], directory=directory
        )
