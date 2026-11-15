"""Gauge, wall-energy and fixed-band wall-amplitude checks on native outputs."""

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

pytest.importorskip("netCDF4")
from generate_ion_energy_restarts import generate
from netCDF4 import Dataset
from scipy.integrate import cumulative_trapezoid

ROOT = Path(__file__).resolve().parents[2]
INPUT = Path(__file__).parent / "data/ion-energy/BOUT.inp"


@pytest.fixture
def hermes_executable():
    executable = Path(
        os.environ.get("HERMES_EXECUTABLE", ROOT / "build/hermes-3")
    ).resolve()
    if not executable.is_file():
        pytest.skip("Build Hermes-3 or set HERMES_EXECUTABLE to run native regressions")
    return executable


def run(executable, folder, *overrides):
    folder.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(INPUT, folder / "BOUT.inp")
    command = [str(executable), "-d", str(folder), *overrides]
    with (folder / "run.log").open("w") as log:
        completed = subprocess.run(
            command,
            cwd=folder,
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=1200,
            check=False,
        )
    assert completed.returncode == 0, (folder / "run.log").read_text()[-6000:]
    return folder / "BOUT.dmp.0.nc"


def gauge_outputs(path):
    with Dataset(path) as dump:
        names = (
            "ddt(Eion)",
            "ddt(Pe)",
            "ddt(Nh+)",
            "Tion",
            "ion_KE_perp",
            "polarisation_flux_z",
            "DivJpol",
            "DivJReynolds",
            "DivJviscous",
            "Pe",
            "ion_pressure_potential",
        )
        fields = {name: np.asarray(dump[name][0, 2:-2, 0, :]) for name in names}
        fields["polarisation_flux_x"] = np.asarray(
            dump["polarisation_flux_x"][0, 1:-2, 0, :]
        )
        nz = fields["Pe"].shape[-1]
        dz = float(np.asarray(dump["dz"][:]).flat[0])
        rho_s = float(dump["rho_s0"][...])
        wave = 2 * np.pi * np.fft.rfftfreq(nz, d=dz)
        pressure = fields["Pe"] + fields["ion_pressure_potential"]
        dia = (
            2
            * rho_s**2
            / 1.5**2
            * np.fft.irfft(1j * wave * np.fft.rfft(pressure, axis=-1), n=nz, axis=-1)
        )
        # The vorticity RHS is the sum of native current diagnostics and the
        # analytic curvature divergence for the uniform identity-coordinate slab.
        fields["ddt(Vort)"] = fields["DivJReynolds"] + fields["DivJviscous"] + dia
        residual = float(np.max(abs(fields["DivJpol"] + dia)))
        assert_allclose(residual, float(dump["ion_chi_residual"][0]), atol=1e-15)
    invariant = (
        "ddt(Eion)",
        "ddt(Pe)",
        "ddt(Nh+)",
        "ddt(Vort)",
        "Tion",
        "ion_KE_perp",
        "polarisation_flux_x",
        "polarisation_flux_z",
    )
    return {key: fields[key] for key in invariant}, residual


@pytest.mark.parametrize("viscosity", (0.0, 1.0))
def test_potential_gauge_invariance(hermes_executable, tmp_path, viscosity):
    outputs = {}
    reference_state = None
    for shift in (-1.0, 0.0, 1.0):
        folder = tmp_path / f"shift{shift:g}"
        restart = generate(folder, shift)
        with Dataset(restart) as dump:
            state = {
                name: np.asarray(dump[name][2:-2, 0, :])
                for name in ("Vort", "Nh+", "Pe", "Eion")
            }
        if reference_state is None:
            reference_state = state
        else:
            for name, value in state.items():
                assert_allclose(value, reference_state[name], atol=0, rtol=0)
        path = run(
            hermes_executable,
            folder,
            "restart=true",
            "solver:nout=0",
            f"evolve_vorticity_ion_energy:viscosity={viscosity}",
            "evolve_vorticity_ion_energy:phi_atol=1e-12",
            "evolve_vorticity_ion_energy:phi_rtol=1e-12",
            "evolve_vorticity_ion_energy:vort_atol=1e-13",
            "evolve_vorticity_ion_energy:phi_max_iterations=1000",
        )
        outputs[shift] = gauge_outputs(path)
    baseline, baseline_residual = outputs[0.0]
    for shift in (-1.0, 1.0):
        shifted, residual = outputs[shift]
        assert max(residual, baseline_residual) < 1e-12
        # A single energy RHS evaluation has gauge forcing C*rchi. The floor
        # also covers current inversion roundoff on this small smooth grid.
        tolerance = 1e-12 + 1.1 * abs(shift) * max(residual, baseline_residual)
        for name, field in shifted.items():
            assert_allclose(field, baseline[name], atol=tolerance, rtol=0, err_msg=name)


def test_total_energy_closes_against_wall_power(hermes_executable, tmp_path):
    path = run(hermes_executable, tmp_path / "wall-energy")
    with Dataset(path) as dump:
        assert int(dump["evolve_vorticity_ion_energy_version"][...]) == 4
        physical = lambda name: np.asarray(dump[name][:, 2:-2, 0, :])
        volume = (
            np.asarray(dump["J"][:])
            * np.asarray(dump["dx"][:])
            * np.asarray(dump["dy"][:])
            * np.asarray(dump["dz"][:])
        )
        if volume.ndim == 2:
            volume = volume[..., None]
        volume = np.broadcast_to(volume, dump["Vort"].shape[1:])[2:-2, 0, :]
        inventory = ((physical("Eion") + 1.5 * physical("Pe")) * volume).sum(
            axis=(1, 2)
        )
        phi = np.asarray(dump["phi"][:, :, 0, :]) + np.asarray(
            dump["ion_pressure_potential"][:, :, 0, :]
        )
        current = np.asarray(dump["polarisation_flux_x"][:, :, 0, :])
        weights = np.asarray(dump["polarisation_wall_energy_weight"][:, :, 0, :])
        work = np.asarray(dump["polarisation_wall_work_weight"][:, :, 0, :])
        wall_power = (
            current[:, -3] * weights[:, -3] - current[:, 1] * weights[:, 1]
        ).sum(axis=-1)
        for face in (1, -3):
            midpoint = (phi[:, face] + phi[:, face + 1]) / 2
            assert_allclose(work[:, face], midpoint, atol=2e-12, rtol=0)
            assert_allclose(
                (current[:, face] * work[:, face]).sum(axis=-1),
                midpoint.mean(axis=-1) * current[:, face].sum(axis=-1),
                atol=1e-15,
                rtol=1e-12,
            )
        terms = physical("DivEpol") * volume
        cancellation_scale = np.max(np.sum(abs(terms), axis=(1, 2)))
        assert_allclose(
            terms.sum(axis=(1, 2)),
            wall_power,
            atol=64 * np.finfo(float).eps * cancellation_scale,
            rtol=1e-10,
        )
        supplied = cumulative_trapezoid(
            wall_power, x=np.asarray(dump["t_array"][:]), initial=0
        )
        assert np.max(abs(inventory - inventory[0] + supplied)) / inventory[0] < 1e-10
        assert physical("Tion").min() > 0
        assert physical("Nh+").min() > 0


@pytest.mark.slow
def test_drift_plane_fixed_band_wall_amplitudes(hermes_executable, tmp_path):
    from energy_analysis.spectra import analyse_native_bands

    folder = tmp_path / "drift-plane"
    run(
        hermes_executable,
        folder,
        "mesh:nx=132",
        "mesh:nz=128",
        "solver:nout=240",
        "solver:output_step=5",
        "solver:atol=1e-10",
        "solver:rtol=1e-5",
        "e:diagnose=false",
        "h+:diagnose=false",
    )
    report = analyse_native_bands(folder)
    assert report["endpoint"] == 1200
    # Rounded bounds include margin above the maximum individual-wall RMS
    # over the complete 128² trajectory for these expressions and tolerances.
    for band, threshold in (("16-31", 3e-5), ("32-63", 5e-6)):
        assert report["vorticity"][band]["maximum_individual_wall_RMS"] < threshold
