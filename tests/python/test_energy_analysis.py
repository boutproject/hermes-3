"""Independent operator identities and inventory/data-contract regression tests."""

import numpy as np
import pytest
import xarray as xr
from boutdata.data import BoutOptions
from energy_analysis import (
    AuditError,
    cadence_comparison,
    compare_diagnostics,
    electron_curvature_work,
    energy_balance,
    energy_densities,
    energy_inventory,
    integrated_power,
    open_run,
    perpendicular_operator,
    polarisation_energy_coefficient,
    regular_timestep,
    spectral_ddz,
)
from numpy.testing import assert_allclose


def fields(nx=16, nz=16, mxg=2):
    x = 2 * np.pi * (np.arange(nx + 2 * mxg) - mxg) / nx
    z = 2 * np.pi * np.arange(nz) / nz
    f = xr.DataArray(
        np.cos(x[:, None] + z[None, :]),
        dims=("x", "z"),
        coords={"x": np.arange(len(x)), "z": np.arange(nz)},
    )
    g = xr.Dataset(
        {
            name: xr.ones_like(f) * value
            for name, value in {
                "J": 1.0,
                "dx": 1.0,
                "dz": 1.0,
                "g11": 1.2,
                "g33": 1.4,
                "g13": 0.3,
            }.items()
        }
    )
    return f, g


def test_fourier_mode_including_cross_metric():
    f, g = fields()
    op = perpendicular_operator(1.5, f, g)
    k = 2 * np.pi / 16
    eigenvalue = -1.5 * (
        4 * (1.2 + 1.4) * np.sin(k / 2) ** 2 + 2 * 0.3 * np.sin(k) ** 2
    )
    assert_allclose(
        op.divergence.isel(x=slice(2, -2)),
        eigenvalue * f.isel(x=slice(2, -2)),
        atol=2e-15,
    )
    assert_allclose(
        op.kperp.isel(x=slice(2, -2)).sum(), -eigenvalue * 16 * 16 / 4, atol=2e-13
    )


def test_constant_and_singleton_z():
    f, g = fields(nz=1)
    op = perpendicular_operator(2, xr.ones_like(f), g)
    assert_allclose(op.divergence.isel(x=slice(2, -2)), 0)
    assert_allclose(op.kperp.isel(x=slice(2, -2)), 0)


def test_nonuniform_metric_symmetry_and_integration_by_parts():
    f, g = fields()
    phase = xr.DataArray(
        2 * np.pi * (np.arange(f.sizes["x"]) - 2) / 16, dims="x", coords={"x": f.x}
    )
    g["dx"] = 1 + 0.2 * np.cos(phase)
    g["dz"] = 1 + 0.1 * np.sin(phase + 2 * np.pi * f.z / 16)
    g["J"] = 1.3 + 0.1 * np.cos(phase)
    g["g13"] = 0.3 + 0.1 * np.cos(phase + 2 * np.pi * f.z / 16)
    a = 1 + 0.2 * np.sin(phase - 2 * np.pi * f.z / 16)
    h = np.sin(phase + 4 * np.pi * f.z / 16)
    lf = perpendicular_operator(a, f, g)
    lh = perpendicular_operator(a, h, g)
    volume = g.J * g.dx * g.dz
    integral = lambda v: float((v * volume).isel(x=slice(2, -2)).sum())
    assert_allclose(integral(f * lf.divergence), -2 * integral(lf.kperp), atol=3e-14)
    assert_allclose(
        integral(h * lf.divergence), integral(f * lh.divergence), atol=3e-14
    )
    assert_allclose(integral(lf.divergence), 0, atol=1e-14)
    # Spatially chunked calculations must preserve neighbor overlaps.
    chunked = perpendicular_operator(
        a.chunk({"x": 5, "z": 7}), f.chunk({"x": 5, "z": 7}), g
    )
    xr.testing.assert_allclose(chunked.compute(), lf)


def test_radial_boundary_work_and_flux_closure():
    f, g = fields()
    f = f + 0.1 * f.x**2
    op = perpendicular_operator(1.0, f, g)
    volume = g.J * g.dx * g.dz
    interior = lambda v: v.isel(x=slice(2, -2))
    work = (
        ((f.isel(x=-3) + f.isel(x=-2)) / 2) * op.flux_x.isel(x=-3)
        - ((f.isel(x=1) + f.isel(x=2)) / 2) * op.flux_x.isel(x=1)
    ).sum()
    assert_allclose(
        (interior((f * op.divergence + 2 * op.kperp) * volume)).sum(), work, atol=3e-12
    )
    outward = (op.flux_x.isel(x=-3) - op.flux_x.isel(x=1)).sum()
    assert_allclose(interior(op.divergence * volume).sum(), outward, atol=2e-13)
    closed = perpendicular_operator(1.0, f, g, boundary_flux=False)
    assert_allclose(interior(closed.divergence * volume).sum(), 0, atol=1e-13)
    assert_allclose(closed.flux_x.isel(x=[1, -3]), 0)


@pytest.fixture
def dataset():
    f, g = fields(nx=8, nz=8)
    ds = g.expand_dims(y=[0]).copy()
    base = (
        xr.ones_like(f)
        .expand_dims(t=[0.0, 2.0, 4.0], y=[0])
        .transpose("t", "x", "y", "z")
    )
    for name, value in {
        "Nd+": 0.5,
        "Nimpurity": 0.25,
        "Ne": 1.0,
        "Tion": 1.0,
        "Pe": 1.0,
        "Eion": 1.125 + 0.006666666666666667,
        "NVd+": 0.1,
        "NVimpurity": -0.05,
        "NVe": 0.01,
        "phi": 0.0,
        "Vort": 0.0,
        "Bxy": 1.0,
        "dy": 1.0,
    }.items():
        ds[name] = (
            base * value if name not in ("Bxy", "dy") else xr.ones_like(ds.J) * value
        )
    options = BoutOptions()
    options["hermes:components"] = "e, d+, impurity, evolve_vorticity_ion_energy"
    for name, mass, charge in (("e", "1/1836", -1), ("d+", 2, 1), ("impurity", 3, 2)):
        options.getSection(name)
        options[name]["AA"] = mass
        options[name]["charge"] = charge
        options[name]["type"] = "evolve_momentum"
        options[name]["density_floor"] = 1e-7
    options["evolve_vorticity_ion_energy:average_atomic_mass"] = "7/4"
    options["evolve_vorticity_ion_energy"]["diagnose"] = "false"
    ds.attrs = {
        "options": options,
        "metadata": {
            "MXG": 2,
            "MYG": 0,
            "nx": 12,
            "ny": 1,
            "keep_xboundaries": True,
            "keep_yboundaries": True,
            "rho_s0": 0.02,
            "Omega_ci": 1e6,
            "Nnorm": 2e18,
            "Tnorm": 5,
        },
    }
    for name in ds.variables:
        if name not in ("x", "y", "z"):
            ds[name].attrs = {"units_type": "hermes", "conversion": 7.0}
    return ds


@pytest.mark.parametrize("settings_present", [False, True])
def test_path_loader_prefers_settings_and_accepts_input_only(
    dataset, tmp_path, monkeypatch, settings_present
):
    import xhermes

    (tmp_path / "BOUT.inp").write_text("supplied input")
    if settings_present:
        (tmp_path / "BOUT.settings").write_text("resolved settings")
    calls = []

    def load(**kwargs):
        calls.append(kwargs)
        return dataset

    monkeypatch.setattr(xhermes, "open_hermesdataset", load)
    run = open_run(tmp_path)
    assert calls[0]["inputfilepath"] == str(
        tmp_path / ("BOUT.settings" if settings_present else "BOUT.inp")
    )
    assert calls[0]["keep_xboundaries"] and calls[0]["keep_yboundaries"]
    assert run.source == str(tmp_path)


def test_multispecies_inventory_units_and_no_double_counting(dataset):
    r = open_run(dataset)
    d = energy_densities(r)
    assert [s.name for s in r.species] == ["e", "d+", "impurity"]
    assert not dataset.options["evolve_vorticity_ion_energy"].get_bool("diagnose")
    assert_allclose(r.interior(d.pressure_potential), (2 * 0.5 + 3 * 0.25 / 2) / 1.75)
    assert_allclose(
        r.interior(d.ion_parallel), 0.1**2 / (2 * 2 * 0.5) + 0.05**2 / (2 * 3 * 0.25)
    )
    assert_allclose(r.interior(d.electron_parallel), 0.01**2 / (2 * (1 / 1836)))
    inventory = energy_inventory(r, densities=d)
    from scipy.constants import elementary_charge

    factor = elementary_charge * 2e18 * 5 * 0.02**3 * 8 * 8
    assert_allclose(inventory.ion_saved, dataset.Eion.isel(x=2, y=0, z=0) * factor)
    assert_allclose(
        inventory.plasma_saved,
        (dataset.Eion.isel(x=2, y=0, z=0) + 1.5 + 0.01**2 * 1836 / 2) * factor,
    )
    assert inventory.ion_saved.attrs["units"] == "J"
    assert_allclose(inventory.t, [0, 2e-6, 4e-6])
    assert not any(v["available"] for v in compare_diagnostics(r).values())


def test_physical_and_normalised_views_match_without_input_mutation(dataset):
    physical = dataset.copy(deep=True)
    for name in physical.variables:
        if name not in ("x", "y", "z", "t"):
            attrs = dict(physical[name].attrs)
            physical[name] = physical[name] * attrs["conversion"]
            physical[name].attrs = {**attrs, "units_type": "SI"}
    attrs = dict(physical.t.attrs)
    physical = physical.assign_coords(t=physical.t * attrs["conversion"])
    physical.t.attrs = {**attrs, "units_type": "SI"}
    before = physical.Eion.copy(deep=True)
    xr.testing.assert_allclose(
        energy_inventory(open_run(physical)), energy_inventory(open_run(dataset))
    )
    xr.testing.assert_identical(before, physical.Eion)


def test_diagnostics_optional_and_cross_checked(dataset):
    r = open_run(dataset)
    d = energy_densities(r)
    for field, reservoir in (
        ("ion_KE_perp", "ion_perpendicular"),
        ("ion_KE_parallel", "ion_parallel"),
        ("ion_pressure_potential", "pressure_potential"),
    ):
        dataset[field] = d[reservoir]
    for value in compare_diagnostics(open_run(dataset)).values():
        assert value["available"] and value["max_absolute_error_normalised"] == 0
    stripped = dataset.drop_vars(
        ["ion_KE_perp", "ion_KE_parallel", "ion_pressure_potential"]
    )
    xr.testing.assert_allclose(
        energy_inventory(open_run(dataset)), energy_inventory(open_run(stripped))
    )


@pytest.mark.parametrize(
    "change,match",
    [
        (lambda ds: ds.drop_vars("Nimpurity"), "Missing density"),
        (lambda ds: ds.drop_vars("g13"), "Missing inventory fields"),
    ],
)
def test_missing_contract(dataset, change, match):
    changed = change(dataset)
    with pytest.raises(AuditError, match=match):
        open_run(changed)


def test_floors_guards_unknown_units_and_momentum(dataset):
    near_floor = dataset.copy(deep=True)
    near_floor["Nd+"][:] = 1e-7
    with pytest.raises(AuditError, match="near floor"):
        open_run(near_floor)
    with pytest.raises(AuditError, match="Missing saved momentum"):
        energy_inventory(open_run(dataset.drop_vars("NVd+")))
    dataset.attrs["metadata"]["keep_xboundaries"] = False
    with pytest.raises(AuditError, match="radial guards"):
        open_run(dataset)
    dataset.attrs["metadata"]["keep_xboundaries"] = True
    dataset.phi.attrs.pop("units_type")
    with pytest.raises(AuditError, match="Unknown normalisation"):
        open_run(dataset)


def test_cumulative_balance_regular_time_and_coverage(dataset):
    inventory = energy_inventory(open_run(dataset))
    q = xr.ones_like(inventory.t) * 2
    q.attrs["units"] = "W"
    zero = q * 0
    zero.attrs["units"] = "W"
    inventory["plasma_saved"] = 3 + 2 * inventory.t
    inventory.plasma_saved.attrs["units"] = "J"
    balance = energy_balance(inventory, source_power=q, outward_power=zero)
    assert balance.attrs["complete"]
    assert_allclose(balance.residual, 0, atol=3e-16)
    assert_allclose(balance.interval_residual.isel(t=slice(1, None)), 0, atol=5e-16)
    partial = energy_balance(inventory)
    assert not partial.attrs["complete"] and "residual" not in partial
    with pytest.raises(AuditError, match="expressed in W"):
        energy_balance(inventory, source_power=q.assign_attrs(units="J"))


@pytest.mark.parametrize("time", [[0], [0, np.nan], [0, 0, 1], [0, 2, 1], [0, 1, 3]])
def test_invalid_time_axes(time):
    with pytest.raises(AuditError):
        regular_timestep(time)


def test_regular_time_roundoff():
    assert_allclose(regular_timestep([1e6, 1e6 + 0.1, 1e6 + 0.2]), 0.1, atol=1e-10)


def test_y_guards_are_excluded(dataset):
    expanded = dataset.isel(y=0, drop=True).expand_dims(y=np.arange(5)).copy(deep=True)
    expanded.attrs["metadata"] = dict(dataset.attrs["metadata"], MYG=2)
    inventory = energy_inventory(open_run(expanded))
    xr.testing.assert_allclose(inventory, energy_inventory(open_run(dataset)))
    # Saved boundary energies must not pollute the physical-cell inventory.
    expanded.Eion.loc[{"y": [0, 1, 3, 4]}] = 1e6
    xr.testing.assert_allclose(energy_inventory(open_run(expanded)), inventory)


def test_common_temperature_floor_is_checked(dataset):
    dataset.options["evolve_vorticity_ion_energy:density_floor"] = 1
    with pytest.raises(AuditError, match="common-temperature recovery floor"):
        open_run(dataset)


def test_algebraic_electron_closure_does_not_add_inertia(dataset):
    dataset.options["e:type"] = "quasineutral, evolve_pressure, zero_current"
    densities = energy_densities(open_run(dataset))
    assert float(densities.electron_parallel.max()) == 0


def test_nonfinite_guard_cannot_silently_pass(dataset):
    dataset.phi.loc[{"x": 1}] = np.nan
    inventory = energy_inventory(open_run(dataset))
    assert bool(np.isnan(inventory.ion_perpendicular).all())


def test_spectral_ddz_modes_nyquist_and_chunks():
    z = 2 * np.pi * np.arange(16) / 16
    field = xr.DataArray(np.sin(3 * z) + 0.2 * (-1.0) ** np.arange(16), dims="z")
    assert_allclose(spectral_ddz(field, 2 * np.pi / 16), 3 * np.cos(3 * z), atol=2e-14)
    xr.testing.assert_allclose(
        spectral_ddz(field.chunk({"z": 5}), 2 * np.pi / 16).compute(),
        spectral_ddz(field, 2 * np.pi / 16),
    )
    assert_allclose(spectral_ddz(field.isel(z=slice(0, 1)), 1), 0)
    with pytest.raises(AuditError, match="uniform positive"):
        spectral_ddz(field, [1, 2])


def test_curvature_work_units_and_diagnostic_independence(dataset):
    options = dataset.options
    for key, value in {"bxcvx": 0, "bxcvy": 0, "bxcvz": 2}.items():
        options["mesh:" + key] = value
    options["mesh:ddz:first"] = "FFT"
    options["mesh:ddz:fft_filter"] = 0
    dataset["dz"] = xr.ones_like(dataset.dz) * (2 * np.pi / 8)
    dataset.dz.attrs["units_type"] = "hermes"
    z = 2 * np.pi * dataset.z / 8
    dataset["phi"] = xr.ones_like(dataset.phi) * np.sin(z)
    dataset.phi.attrs["units_type"] = "hermes"
    dataset["Pe"] = xr.ones_like(dataset.Pe) * np.cos(z)
    dataset.Pe.attrs["units_type"] = "hermes"
    run = open_run(dataset)
    work = electron_curvature_work(run)
    assert_allclose(
        work, -4 * run.length_scale**2 * xr.ones_like(work) * np.cos(z) ** 2, atol=1e-18
    )
    power = integrated_power(run, work)
    assert power.attrs["units"] == "W"
    assert_allclose(
        power,
        run.integrate(work) * run.energy_scale * run.length_scale**3 * run.frequency,
    )
    assert_allclose(power.t, dataset.t / run.frequency)
    options["mesh:ddz:fft_filter"] = 0.1
    with pytest.raises(AuditError, match="unfiltered"):
        electron_curvature_work(run)
    options["mesh:ddz:fft_filter"] = 0
    options["mesh:bxcvx"] = 1
    with pytest.raises(AuditError, match="Z-only"):
        electron_curvature_work(run)


def test_multispecies_polarisation_energy_coefficient(dataset):
    run = open_run(dataset)
    result = polarisation_energy_coefficient(run)
    rho = 2 * 0.5 + 3 * 0.25
    kinetic = (2 * 0.1**2 / (2 * 2 * 0.5) + (3 / 2) * 0.05**2 / (2 * 3 * 0.25)) / rho
    assert_allclose(result.parallel_kinetic, kinetic)
    assert_allclose(result.enthalpy, 2.5 * (2 * 0.5 + 3 / 2 * 0.25) / rho)
    assert_allclose(result.total, result.enthalpy + result.parallel_kinetic)
    # Actual mass density, rather than the Boussinesq reference, controls weights.
    dataset["Nimpurity"] *= 1.5
    changed = polarisation_energy_coefficient(open_run(dataset))
    assert float(changed.enthalpy.max()) != float(result.enthalpy.max())


def test_cadence_refinement_distinguishes_quadrature_error():
    def budget(dt):
        time = np.arange(0, 1 + dt / 2, dt)
        inventory = xr.Dataset(
            {"electron_thermal": ("t", 2 + time**3 / 3)}, coords={"t": time}
        )
        inventory.electron_thermal.attrs["units"] = "J"
        power = xr.DataArray(
            time**2, dims="t", coords={"t": time}, attrs={"units": "W"}
        )
        zero = xr.zeros_like(power).assign_attrs(units="W")
        return energy_balance(
            inventory,
            reservoir="electron_thermal",
            source_power=power,
            outward_power=zero,
        )

    coarse, fine = budget(0.25), budget(0.125)
    comparison = cadence_comparison(coarse, fine)
    assert comparison["refinement_ratio"] == 2
    assert_allclose(
        comparison["fine_max_residual_J"], comparison["coarse_max_residual_J"] / 4
    )
    assert_allclose(
        comparison["estimated_fine_sampling_error_J"], comparison["fine_max_residual_J"]
    )
    assert comparison["estimated_cadence_independent_max_residual_J"] < 1e-15
    with pytest.raises(AuditError, match="nested saved times"):
        cadence_comparison(coarse, fine.isel(t=slice(1, None)))
    with pytest.raises(AuditError, match="integer refinement"):
        cadence_comparison(coarse, coarse)


def test_transport_flux_gauss_and_weight_linearity():
    from energy_analysis.transport import arakawa, face_divergence, xppm

    f, g = fields()
    g["dy"] = xr.ones_like(g.J)
    stream = np.sin(2 * np.pi * f.x / 16) + 0.3 * np.cos(4 * np.pi * f.z / 16)
    op = xppm(1 + f, stream, g)
    volume = g.J * g.dx * g.dy * g.dz
    assert_allclose((op.divergence * volume).isel(x=slice(2, -2)).sum(), 0, atol=1e-13)
    w = 0.4 + 0.1 * f
    weighted = xppm(1 + f, stream, g, weight=w)
    other = xppm(1 + f, stream, g, weight=1 - w)
    xr.testing.assert_allclose(weighted + other, op)
    xr.testing.assert_allclose(
        face_divergence({"x": op.flux_x, "z": op.flux_z}, g), op.divergence
    )
    h = np.sin(2 * np.pi * (f.x - 2) / 16 - 4 * np.pi * f.z / 16)
    bracket = arakawa(f, h, g.dx, g.dz).isel(x=slice(2, -2))
    assert_allclose(bracket.sum(), 0, atol=1e-14)
    assert_allclose((f.isel(x=slice(2, -2)) * bracket).sum(), 0, atol=1e-14)


def test_parallel_face_flux_and_open_midpoints():
    from energy_analysis.transport import parallel_flux

    y = np.arange(20)
    f = xr.DataArray(
        1 + 0.1 * np.cos(2 * np.pi * (y - 2) / 16), dims="y", coords={"y": y}
    )
    v = 0.3 + 0.05 * np.sin(2 * np.pi * (y - 2) / 16)
    v = xr.DataArray(v, dims="y", coords={"y": y})
    g = xr.Dataset({n: xr.ones_like(f) for n in ("J", "dx", "dy", "dz", "g_22")})
    op = parallel_flux(f, v, xr.ones_like(f), g)
    assert_allclose(op.divergence.isel(y=slice(2, -2)).sum(), 0, atol=1e-15)
    opened = parallel_flux(f, v, xr.ones_like(f), g, periodic_y=False)
    expected = ((f + f.shift(y=-1)) / 2) * ((v + v.shift(y=-1)) / 2)
    assert_allclose(opened.flux_y.isel(y=[1, -3]), expected.isel(y=[1, -3]))
    assert_allclose(
        opened.divergence.isel(y=slice(2, -2)).sum(),
        opened.flux_y.isel(y=-3) - opened.flux_y.isel(y=1),
    )


def test_sparse_polarisation_solve_and_neumann_guards(dataset):
    from energy_analysis.fluxes import fill_transport_guards, solve_polarisation

    dataset["g13"] *= 0
    dataset["g_22"] = xr.ones_like(dataset.J)
    run = open_run(dataset)
    phi = xr.ones_like(dataset.phi) * np.cos(
        2 * np.pi * (dataset.x - 2) / 8 + 2 * np.pi * dataset.z / 8
    )
    rhs = perpendicular_operator(1.75, phi, dataset).divergence
    actual = solve_polarisation(run, rhs, True, True)
    lazy = solve_polarisation(run, rhs.chunk({"t": 1, "x": 3, "z": 3}), True, True)
    xr.testing.assert_allclose(lazy.compute(), actual)
    xr.testing.assert_allclose(run.interior(actual), run.interior(phi), atol=2e-14)
    guarded = fill_transport_guards(phi, run, False, True)
    xr.testing.assert_allclose(guarded.isel(x=1, drop=True), phi.isel(x=2, drop=True))
    xr.testing.assert_allclose(guarded.isel(x=0, drop=True), phi.isel(x=3, drop=True))
    rhs = rhs + 0.1
    with pytest.raises(AuditError, match="incompatible"):
        solve_polarisation(run, rhs, True, True)


def test_prescribed_source_conversion_and_coverage(dataset):
    from energy_analysis.fluxes import PrescribedSources, reconstruct_sources
    from scipy.constants import proton_mass

    run = open_run(dataset)
    zero = {s.name: 0 for s in run.species}
    rates = PrescribedSources(
        {**zero, "d+": 20},
        {**zero, "d+": 3e20},
        {**zero, "d+": 0.1},
        "analytic SI source",
    )
    result = reconstruct_sources(run, rates)
    v = dataset["NVd+"] / (2 * dataset["Nd+"])
    expect = (
        20 / (run.energy_scale * run.frequency)
        + v * 0.1 / (proton_mass * 2e18 * run.length_scale * run.frequency**2)
        - v**2 * 3e20 / (2e18 * run.frequency)
    )
    xr.testing.assert_allclose(result.Eion_src, expect)
    with pytest.raises(AuditError, match="every species"):
        reconstruct_sources(run, PrescribedSources({}, zero, zero, "incomplete"))


def test_parallel_numerical_heating_kinetic_flux_identity():
    from energy_analysis.transport import (
        numerical_heating_defect,
        numerical_parallel_heating,
        parallel_flux,
        parallel_kinetic_flux,
    )

    y = np.arange(36)
    angle = 2 * np.pi * (y - 2) / 32
    n = xr.DataArray(1 + 0.3 * np.cos(angle), dims="y", coords={"y": y})
    v = xr.DataArray(0.4 + 0.2 * np.sin(2 * angle), dims="y", coords={"y": y})
    wave = xr.ones_like(n) * 0.8
    geometry = xr.Dataset(
        {key: xr.ones_like(n) for key in ("J", "dx", "dy", "dz", "g_22")}
    )
    dn = -parallel_flux(n, v, wave, geometry).divergence
    dnv = -parallel_flux(n, v, wave, geometry, momentum=True).divergence
    kinetic = v * dnv - 0.5 * v**2 * dn
    heating = numerical_parallel_heating(n, v, wave, geometry)
    divergence = parallel_kinetic_flux(n, v, wave, geometry).divergence
    assert_allclose(
        (kinetic + heating + divergence).isel(y=slice(2, -2)),
        numerical_heating_defect(n, v, wave, geometry).isel(y=slice(2, -2)),
        atol=2e-16,
    )
    assert (
        float(numerical_heating_defect(n, v, wave, geometry).isel(y=slice(2, -2)).sum())
        > 1e-4
    )


def test_flux_report_explicitly_incomplete_for_unsupported_metrics(dataset):
    from energy_analysis.report import flux_report

    report = flux_report(open_run(dataset))
    assert not report["flux_audit"]["complete"]
    assert "requires saved g_22" in report["flux_audit"]["reason"]


def test_open_rectangle_is_not_misclassified_as_tokamak(dataset):
    metadata = dataset.attrs["metadata"]
    metadata.update(ixseps1=0, ixseps2=0, jyseps1_1=-1, jyseps2_2=0)
    run = open_run(dataset)
    assert run.myg == 0
    metadata["ixseps1"] = 4
    with pytest.raises(AuditError, match="Tokamak boundary masks"):
        open_run(dataset)


def test_sparse_dirichlet_solve_before_transport_boundary_change(dataset):
    from energy_analysis.fluxes import solve_polarisation

    dataset["g13"] *= 0
    run = open_run(dataset)
    nx = dataset.sizes["x"] - 4
    xphase = np.pi * (dataset.x - 2 + 0.5) / nx
    phi = xr.ones_like(dataset.phi) * np.sin(xphase) * np.cos(2 * np.pi * dataset.z / 8)
    rhs = perpendicular_operator(1.75, phi, dataset).divergence
    solved = solve_polarisation(run, rhs, False, True)
    xr.testing.assert_allclose(run.interior(solved), run.interior(phi), atol=2e-14)
    # The transport potential is Neumann, although the inverse used Dirichlet.
    xr.testing.assert_allclose(solved.isel(x=1, drop=True), solved.isel(x=2, drop=True))


def test_relaxation_auxiliary_solve_retains_full_current_boundary(dataset):
    from energy_analysis.fluxes import solve_polarisation
    from energy_analysis.transport import fv_xz

    dataset["g13"] *= 0
    dataset.options["evolve_vorticity_ion_energy:phi_boundary_relax"] = "true"
    run = open_run(dataset)
    nx = dataset.sizes["x"] - 4
    potential = (
        xr.ones_like(dataset.phi)
        * np.sin(np.pi * (dataset.x - 2 + 0.5) / nx)
        * np.cos(2 * np.pi * dataset.z / 8)
    )
    rhs = perpendicular_operator(1.75, potential, dataset).divergence
    solved = solve_polarisation(run, rhs, False, True)
    xr.testing.assert_allclose(
        run.interior(solved), run.interior(potential), atol=2e-14
    )
    xr.testing.assert_allclose(
        solved.isel(x=1, drop=True), -solved.isel(x=2, drop=True)
    )
    full = fv_xz(1.75, solved, dataset)
    xr.testing.assert_allclose(
        run.interior(full.divergence), run.interior(rhs), atol=2e-14
    )
    assert float(abs(full.flux_x.isel(x=[1, -3])).max()) > 1e-3
    for weight in (1.0, 0.7, -0.2, 2.5, 0.03, 0.4):
        closed = fv_xz(weight, solved, dataset, boundary_flux=False)
        assert_allclose(closed.flux_x.isel(x=[1, -3]), 0)


def test_boundary_relaxation_core_mean_excludes_y_guards(dataset):
    from energy_analysis.fluxes import boundary_relaxation

    ds = dataset.isel(y=0, drop=True).expand_dims(y=np.arange(8)).copy()
    for name in ("g13", "g23", "g_23"):
        ds[name] = xr.zeros_like(ds.J)
    ds["g_22"] = xr.ones_like(ds.J)
    opts = ds.options
    opts["mesh:paralleltransform:type"] = "identity"
    opts["mesh:ddz:first"] = "FFT"
    opts["mesh:ddz:fft_filter"] = 0
    opts["mesh:bxcvz"] = 0
    opts["periodicX"] = "false"
    opts["evolve_vorticity_ion_energy:poloidal_flows"] = "false"
    opts["evolve_vorticity_ion_energy:phi_boundary_relax"] = "true"
    opts["evolve_vorticity_ion_energy:phi_core_averagey"] = "true"
    opts["evolve_vorticity_ion_energy:phi_boundary_timescale"] = "2e-6"
    for name in ("e", "d+", "impurity"):
        opts[name]["use_div_par_fvv"] = "true"
    ds.attrs["metadata"].update(
        MYG=2,
        keep_yboundaries=True,
        nx=12,
        ny=4,
        ixseps1=12,
        jyseps1_1=-1,
        jyseps2_2=3,
        HERMES_SLOPE_LIMITER="MC",
    )
    interior = xr.where((ds.y >= 2) & (ds.y < 6), 1 + 0.1 * ds.y, 999)
    ds["phi"] = xr.ones_like(ds.phi) * interior
    for guard, cell in ((1, 2), (-2, -3)):
        ds.phi[{"x": guard}] = 0.4 - ds.phi.isel(x=cell)
    ds.phi.attrs = dataset.phi.attrs
    run = open_run(ds)
    boundary = boundary_relaxation(run)
    assert_allclose(boundary.phi_boundary, 0.2)
    assert_allclose(boundary.phi_boundary_target.sel(boundary="inner"), 1.35)
    assert_allclose(boundary.phi_boundary_rate.sel(boundary="inner"), (1.35 - 0.2) / 2)
    assert_allclose(
        boundary.phi_boundary_target.sel(boundary="outer"),
        xr.ones_like(boundary.phi_boundary_target.sel(boundary="outer"))
        * (1 + 0.1 * boundary.y),
    )


def test_flux_velocity_requires_saved_electron_closure(dataset):
    from energy_analysis.fluxes import species_velocity

    dataset.options["e:type"] = "quasineutral, evolve_pressure, zero_current"
    run = open_run(dataset.drop_vars("NVe"))
    electron = next(s for s in run.species if s.name == "e")
    with pytest.raises(AuditError, match="flowing species e"):
        species_velocity(run, electron)


def test_sparse_solve_handles_auxiliary_geometry_coordinates(dataset):
    from energy_analysis.fluxes import solve_polarisation

    dataset["g13"] *= 0
    dataset = dataset.set_coords("dx")
    run = open_run(dataset)
    potential = xr.ones_like(dataset.phi) * np.cos(
        2 * np.pi * (dataset.x - 2) / 8 + 2 * np.pi * dataset.z / 8
    )
    rhs = perpendicular_operator(1.75, potential, dataset).divergence
    solved = solve_polarisation(run, rhs, True, True)
    # Reindexing the physical solve must not create conflicting NaN dx guards.
    merged = xr.Dataset({"solved": solved, "original": potential})
    xr.testing.assert_allclose(merged.dx, dataset.dx)


def test_current_face_flux_is_not_labelled_as_energy_power(dataset):
    from energy_analysis import flux_power

    run = open_run(dataset)
    flux = xr.Dataset(
        {"DivJpol": xr.zeros_like(dataset.phi).assign_attrs(quantity="current")}
    )
    with pytest.raises(AuditError, match="energy-flux family"):
        flux_power(run, flux, "DivJpol")


def test_pressure_adapter_has_native_and_common_but_no_saved_ion_energy(dataset):
    ds = dataset.drop_vars(["Eion", "Tion"]).copy(deep=True)
    ds["g13"] = xr.zeros_like(ds.g13)
    ds.g13.attrs = {"units_type": "hermes"}
    options = ds.options
    options["hermes:components"] = "e, d+, impurity, vorticity, polarisation_drift"
    options["vorticity:average_atomic_mass"] = "7/4"
    options["vorticity:diamagnetic_polarisation"] = "true"
    for name in ("d+", "impurity"):
        ds["P" + name] = ds["N" + name].copy()
    ds["Pe"] = ds.Pe * (1 + 0.1 * np.sin(2 * np.pi * ds.z / 8))
    ds.Pe.attrs = {"units_type": "hermes"}
    run = open_run(ds)
    densities = energy_densities(run)
    assert run.formulation == "pressure"
    assert "ion_saved" not in densities
    assert "plasma_saved" not in densities
    # Electron gradients change the native pressure potential, but not the
    # common ion-only FV estimator.
    assert float(abs(run.interior(densities.ion_perpendicular)).max()) < 1e-14
    assert float(abs(run.interior(densities.ion_perpendicular_native)).max()) > 0
    assert "plasma_native" in energy_inventory(run)
