"""Native radial-wall convergence and local multispecies transport identities."""

from pathlib import Path

import numpy as np
from energy_analysis import open_run
from energy_analysis.operators import face_divergence, weighted_faces
from energy_analysis.transfers import (
    polarisation_energy_coefficient,
    polarisation_wall_weights,
)
from numpy.testing import assert_allclose
from test_ion_energy_native import run

DATA = Path(__file__).parent / "data"


def check_wall_work(case):
    ds = case.data
    power = np.zeros(ds.sizes["t"])
    for face, sign in ((case.mxg - 1, -1), (ds.sizes["x"] - case.mxg - 1, 1)):
        selected = {"x": face, "y": slice(case.myg, -case.myg or None)}
        current = ds.polarisation_flux_x.isel(selected)
        work = ds.polarisation_wall_work_weight.isel(selected)
        midpoint = (
            (ds.phi + ds.ion_pressure_potential).isel(x=face, drop=True)
            + (ds.phi + ds.ion_pressure_potential).isel(x=face + 1, drop=True)
        ) / 2
        midpoint = midpoint.isel(y=selected["y"])
        assert_allclose(current * (work - midpoint), 0, atol=1e-18, rtol=0)
        power += sign * np.asarray(
            (current * ds.polarisation_wall_energy_weight.isel(selected)).sum(
                ("y", "z")
            )
        )
    rate = case.integrate(ds["ddt(Eion)"] + 1.5 * ds["ddt(Pe)"])
    scale = case.integrate(abs(ds["ddt(Eion)"]) + 1.5 * abs(ds["ddt(Pe)"])) + abs(power)
    assert float((abs(rate + power) / scale.clip(min=1e-300)).max()) < 1e-11
    assert float(case.interior(ds.Tion).min()) > 0


def test_exact_frozen_wall_rates_converge(hermes_executable, tmp_path):
    rates = []
    for n in (16, 32, 64):
        folder = tmp_path / str(n)
        run(
            hermes_executable,
            folder,
            f"mesh:nx={n + 4}",
            f"mesh:nz={n}",
            input_path=DATA / "ion-energy-frozen/BOUT.inp",
        )
        case = open_run(folder)
        try:
            ds = case.data
            assert int(ds.metadata["NXPE"]) == 1
            x = (np.arange(n) + 0.5) / n
            z = 2 * np.pi * np.arange(n) / n
            pressure = 1 + 0.05 * np.cos(np.pi * x)[:, None] * np.sin(z)
            phi = (
                0.01 * np.sin(np.pi * x)[:, None] * np.cos(z)
                + 0.01 * np.sin(2 * np.pi * x)[:, None]
            )
            physical = lambda f, case=case: np.asarray(
                case.interior(f).isel(t=0, y=0).transpose("x", "z")
            )
            assert_allclose(
                physical(ds.ion_pressure_potential), pressure, atol=1e-11, rtol=0
            )
            assert_allclose(physical(ds.phi), phi, atol=1e-11, rtol=0)
            assert float(ds.ion_phi_iterations.max()) <= 200
            wall = ds.isel(x=[case.mxg, ds.sizes["x"] - case.mxg - 1], t=0, y=0)
            assert float(wall.Tion.std("z").min()) > 0.03
            values = np.array(
                [
                    np.sqrt(np.mean(np.asarray(wall[name]) ** 2, axis=-1))
                    for name in ("ddt(Eion)", "DivJReynolds")
                ]
            )
            assert np.all(values > 1e-5)
            rates.append(values)
            check_wall_work(case)
        finally:
            case.data.close()
    differences = abs(np.diff(rates, axis=0))
    factor = differences[0] / differences[1]
    # Both walls and both native rates approach their limits at second order.
    assert np.all((factor > 3.5) & (factor < 4.5)), (rates, differences, factor)


def test_multispecies_local_face_sources(hermes_executable, tmp_path):
    folder = tmp_path / "mixture"
    run(hermes_executable, folder, input_path=DATA / "ion-energy-multispecies/BOUT.inp")
    case = open_run(folder)
    try:
        ds = case.data
        assert_allclose(float(ds.t[-1]), 0.1, atol=1e-14, rtol=0)
        weights = polarisation_wall_weights(case)
        faces = {"x": ds.polarisation_flux_x, "z": ds.polarisation_flux_z}
        rho = 2 * ds["Nd+"] + 4 * ds["Nhe++"]
        for name, mass, charge in (("d+", 2, 1), ("he++", 4, 2)):
            for quantity, field, source in (
                ("particle", "N", "SN"),
                ("momentum", "NV", "SNV"),
            ):
                coefficient = mass * ds[field + name] / (charge * rho)
                flux = weighted_faces(
                    faces,
                    coefficient,
                    mxg=case.mxg,
                    wall_weight=weights[quantity + "_" + name],
                )
                error = case.interior(ds[source + name] + face_divergence(flux, ds))
                assert_allclose(error, 0, atol=1e-12, rtol=0, err_msg=quantity + name)
        energy = weighted_faces(
            faces,
            polarisation_energy_coefficient(case).total,
            mxg=case.mxg,
            wall_weight=weights.total,
        )
        assert_allclose(
            case.interior(ds.DivEpol - face_divergence(energy, ds)),
            0,
            atol=1e-12,
            rtol=0,
        )
        charge_source = ds["SNd+"] + 2 * ds["SNhe++"]
        assert_allclose(
            case.interior(charge_source + ds.DivJpol), 0, atol=1e-12, rtol=0
        )
        assert_allclose(
            case.interior(
                charge_source
                - (ds.ion_current_rhs - ds.DivJReynolds)
                + ds.ion_auxiliary_residual
            ),
            0,
            atol=1e-12,
            rtol=0,
        )
        for face in (case.mxg - 1, ds.sizes["x"] - case.mxg - 1):
            selected = {"x": face, "y": slice(case.myg, -case.myg)}
            fraction = (2 * ds["Nd+"] / rho).isel(selected).isel(t=0)
            assert float(fraction.std("z").min()) > 0.07
            assert_allclose(
                ds.polarisation_wall_energy_weight.isel(selected),
                weights.total.isel(selected),
                atol=1e-12,
                rtol=0,
            )
            assert_allclose(
                ds.polarisation_wall_thermal_weight.isel(selected),
                (weights.thermal + weights.parallel_kinetic).isel(selected),
                atol=1e-12,
                rtol=0,
            )
        check_wall_work(case)
    finally:
        case.data.close()
