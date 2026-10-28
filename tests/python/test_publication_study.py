"""Independent analytic checks for publication estimators and field mapping."""

import importlib.util
from pathlib import Path

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
    source = (
        Path(__file__).resolve().parents[2]
        / "tests/integrated/2D-energy/publication_scan.py"
    )
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
    args = runner.parser().parse_args(["--end-time", "120", "--outputs", "10"])
    assert runner.configuration(args)["outputs"] == 20


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
    with pytest.raises(AuditError, match="Run failed"):
        runner.execute(exe, directory, {"nproc": 1}, 0.1)
    status = json.loads((directory / "status.json").read_text())
    assert status["exit_code"] == -1
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
