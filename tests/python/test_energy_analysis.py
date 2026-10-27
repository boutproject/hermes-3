"""Independent operator identities and inventory/data-contract regression tests."""

import numpy as np
import pytest
import xarray as xr
from boutdata.data import BoutOptions
from energy_analysis import (
    AuditError,
    compare_diagnostics,
    energy_balance,
    energy_densities,
    energy_inventory,
    open_run,
    perpendicular_operator,
    regular_timestep,
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
    balance = energy_balance(inventory, source_power=q, outward_power=zero)
    assert balance.attrs["complete"]
    assert_allclose(balance.residual, 0, atol=3e-16)
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
