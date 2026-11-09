"""Skew Reynolds transport against an independently differentiated covariant form."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_allclose

s = pytest.importorskip("sympy")
from energy_analysis.operators import perpendicular_operator
from energy_analysis.transport import (
    arakawa,
    face_divergence,
    face_work,
    fv_xz,
    reynolds_face_fluxes,
    solenoidal_reynolds_faces,
    xppm_faces,
)

spec = importlib.util.spec_from_file_location(
    "reynolds_reference",
    Path(__file__).parents[1] / "unit/generate_reynolds_reference.py",
)
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)


def sample(case, n):
    h = (1 if case == "annulus" else 2 * np.pi) / n
    xx = (np.arange(n + 4) - 2 + 0.5) * h + (1 if case == "annulus" else 0)
    zz = (np.arange(n) + 0.5) * 2 * np.pi / n
    shape = (n + 4, n)
    values = {}
    for name, expression in reference.reference(case).items():
        value = s.lambdify((reference.x, reference.z), expression, "numpy", cse=True)(
            xx[:, None], zz[None, :]
        )
        values[name] = xr.DataArray(
            np.broadcast_to(value, shape).copy(), dims=("x", "z")
        )
    geometry = xr.Dataset({name: values[name] for name in ("J", "g11", "g13", "g33")})
    geometry["dx"], geometry["dy"], geometry["dz"] = h, 1.0, 2 * np.pi / n
    geometry["Bxy"] = np.sqrt(values["D"])
    return values, geometry


def periodic_halos(field):
    field = field.copy()
    for i in range(2):
        field[{"x": i}] = field.isel(x=field.sizes["x"] - 4 + i)
        field[{"x": field.sizes["x"] - 2 + i}] = field.isel(x=2 + i)
    return field


@pytest.mark.parametrize("case", ("periodic", "annulus"))
@pytest.mark.parametrize("limiter, minimum", ((None, 1.8), ("mc", 1.5)))
def test_covariant_skew_convergence(case, limiter, minimum):
    errors = []
    for n in (32, 64, 128):
        v, g = sample(case, n)
        faces = reynolds_face_fluxes(
            v["c"],
            v["omega"],
            v["phi"],
            v["psi"],
            g,
            walls=case == "annulus",
            periodic_x=case != "annulus",
            limiter=limiter,
            fill_guards=periodic_halos if case != "annulus" else None,
        )
        error = face_divergence(faces, g) - v["expected"]
        interior = error.isel(x=slice(5, -5))
        errors.append(float(np.sqrt((interior**2).mean())))
        boundary = error.isel(x=[2, 3, -4, -3])
        print(
            case,
            limiter,
            n,
            "interior",
            errors[-1],
            "boundary",
            float(np.sqrt((boundary**2).mean())),
            "metric mismatch",
            float(abs((g.g11 * g.g33 - g.g13**2) / g.Bxy**2 - 1).max()),
        )
    orders = np.log2(np.asarray(errors[:-1]) / errors[1:])
    assert np.all(orders >= minimum), (case, limiter, errors, orders)


@pytest.mark.parametrize("case", ("periodic", "annulus"))
@pytest.mark.parametrize("limiter, minimum", ((None, 1.8), ("mc", 1.5)))
def test_covariant_current_face_convergence(case, limiter, minimum):
    # Divergence alone cannot distinguish the signs of a solenoidal curl.
    # Compare J*Gamma at the actual faces with the independent covariant form.
    errors = {dim: [] for dim in ("x", "z")}
    for n in (32, 64, 128):
        v, g = sample(case, n)
        faces = reynolds_face_fluxes(
            v["c"],
            v["omega"],
            v["phi"],
            v["psi"],
            g,
            walls=case == "annulus",
            periodic_x=case != "annulus",
            limiter=limiter,
            fill_guards=periodic_halos if case != "annulus" else None,
        )
        dx, dz = float(g.dx), float(g.dz)
        xx = (np.arange(n + 4) - 2 + 0.5) * dx + (1 if case == "annulus" else 0)
        zz = (np.arange(n) + 0.5) * dz
        for dim in ("x", "z"):
            expression = reference.reference(case)["expected_" + dim]
            expected = s.lambdify(
                (reference.x, reference.z), expression, "numpy", cse=True
            )(
                (xx + (dx / 2 if dim == "x" else 0))[:, None],
                (zz + (dz / 2 if dim == "z" else 0))[None, :],
            )
            area = dz if dim == "x" else dx
            error = (np.asarray(faces[dim]) / area - expected)[5:-5]
            errors[dim].append(float(np.sqrt(np.mean(error**2))))
    for dim, values in errors.items():
        orders = np.log2(np.asarray(values[:-1]) / values[1:])
        print(case, limiter, dim, "face errors", values, "orders", orders)
        assert np.all(orders >= minimum), (case, limiter, dim, values, orders)


def test_curl_work_and_wall_identities():
    v, g = sample("annulus", 32)
    phi, psi = v["phi"].copy(), v["psi"].copy()
    for cell, step in ((2, -1), (phi.sizes["x"] - 3, 1)):
        for layer in (1, 2):
            phi[{"x": cell + step * layer}] = -phi.isel(x=cell)
            psi[{"x": cell + step * layer}] = -psi.isel(x=cell)
    Phi = phi + psi
    curl = solenoidal_reynolds_faces(v["c"], phi, Phi, g)
    assert_allclose(face_divergence(curl, g).isel(x=slice(2, -2)), 0, atol=1e-13)
    faces = reynolds_face_fluxes(v["c"], v["omega"], phi, psi, g)
    for face in (1, phi.sizes["x"] - 3):
        assert_allclose(faces["x"].isel(x=face), curl["x"].isel(x=face), atol=1e-13)
        assert abs(float(faces["x"].isel(x=face).sum())) < 1e-13
    volume = g.J * g.dx * g.dy * g.dz
    integral = lambda f: float((f * volume).isel(x=slice(2, -2)).sum())
    assert_allclose(
        integral(face_work(faces, Phi, g)),
        -integral(Phi * face_divergence(faces, g)),
        rtol=1e-13,
        atol=1e-13,
    )


def test_flat_slab_agrees_with_legacy_split_interior():
    errors = []
    for n in (32, 64, 128):
        v, g = sample("periodic", n)
        g["J"], g["g11"], g["g33"], g["g13"] = 1.0, 1.0, 1.0, 0.0
        phi, psi = v["phi"], v["psi"]
        c = xr.ones_like(phi)
        omega = periodic_halos(perpendicular_operator(c, phi + psi, g).divergence)
        Q = periodic_halos(arakawa(phi, psi, g.dx, g.dz))
        # The old Delp2 uses FD in X and FFT in Z; the new Lambda uses FV in both.
        lapz = s.diff(reference.reference("periodic")["phi"], reference.z, 2)
        h = float(g.dx)
        xx = (np.arange(n + 4) - 2 + 0.5) * h
        zz = (np.arange(n) + 0.5) * h
        analytic_z = xr.DataArray(
            s.lambdify((reference.x, reference.z), lapz, "numpy")(
                xx[:, None], zz[None, :]
            ),
            dims=("x", "z"),
        )
        from energy_analysis.transport import neighbor

        lam = periodic_halos(
            (
                (neighbor(phi, "x") - 2 * phi + neighbor(phi, "x", -1)) / h**2
                + analytic_z
            )
            / 2
        )
        old_a = xppm_faces(omega / 2, phi, g)
        old_b = xppm_faces(lam, phi + psi, g)
        old_q = fv_xz(c / 2, Q, g, periodic_x=True)
        old = {d: old_a[d] + old_b[d] + old_q["flux_" + d] for d in ("x", "z")}
        new = reynolds_face_fluxes(
            c,
            omega,
            phi,
            psi,
            g,
            walls=False,
            periodic_x=True,
            fill_guards=periodic_halos,
        )
        err = face_divergence(new, g) - face_divergence(old, g)
        errors.append(float(np.sqrt((err.isel(x=slice(5, -5)) ** 2).mean())))
    assert np.all(np.log2(np.asarray(errors[:-1]) / errors[1:]) > 1.8), errors
