"""Skew Reynolds transport against an independently differentiated covariant form."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_allclose

s = pytest.importorskip("sympy")
from energy_analysis.fluxes import fill_wall_guards
from energy_analysis.operators import perpendicular_operator
from energy_analysis.transport import (
    face_divergence,
    face_work,
    fv_xz,
    reynolds_face_fluxes,
    weighted_faces,
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


@pytest.mark.parametrize("case", ("periodic", "annulus", "slab"))
def test_covariant_skew_convergence(case):
    minimum = 1.5
    errors = []
    for n in (32, 64, 128):
        v, g = sample(case, n)
        faces = reynolds_face_fluxes(
            v["c"],
            v["omega"],
            v["phi"],
            v["psi"],
            g,
            periodic_x=case != "annulus",
        )
        error = face_divergence(faces, g) - v["expected"]
        interior = error.isel(x=slice(5, -5))
        errors.append(float(abs(interior).mean()))
        boundary = error.isel(x=[2, 3, -4, -3])
        print(
            case,
            "MC vorticity",
            n,
            "interior",
            errors[-1],
            "boundary",
            float(np.sqrt((boundary**2).mean())),
            "metric mismatch",
            float(abs((g.g11 * g.g33 - g.g13**2) / g.Bxy**2 - 1).max()),
        )
    orders = np.log2(np.asarray(errors[:-1]) / errors[1:])
    assert np.all(orders >= minimum), (case, errors, orders)


@pytest.mark.parametrize("case", ("periodic", "annulus", "slab"))
def test_covariant_current_face_convergence(case):
    minimum = 1.5
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
            periodic_x=case != "annulus",
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
            errors[dim].append(float(np.mean(abs(error))))
    for dim, values in errors.items():
        orders = np.log2(np.asarray(values[:-1]) / values[1:])
        print(case, "MC vorticity", dim, "face errors", values, "orders", orders)
        assert np.all(orders >= minimum), (case, dim, values, orders)


def test_reynolds_gauss_and_wall_work():
    v, g = sample("annulus", 32)
    wall = SimpleNamespace(mxg=2)
    phi = fill_wall_guards(v["phi"], wall, "zero")
    psi = fill_wall_guards(v["psi"], wall, "pressure")
    Phi = phi + psi
    faces = reynolds_face_fluxes(v["c"], v["omega"], phi, psi, g)
    volume = g.J * g.dx * g.dy * g.dz
    scale = max(
        1.0,
        float((abs(faces["x"]) + abs(faces["z"])).isel(x=slice(2, -2)).max())
        / float(volume.min()),
    )
    rng = np.random.default_rng(1729)
    weight = xr.DataArray(rng.uniform(0.1, 2, phi.shape), dims=phi.dims)
    changed = weight.copy()
    changed[{"x": [0, 1, -2, -1]}] = 1e6
    first = weighted_faces(faces, weight)
    second = weighted_faces(faces, changed)
    assert_allclose(
        face_divergence(first, g).isel(x=slice(2, -2)),
        face_divergence(second, g).isel(x=slice(2, -2)),
        atol=0,
        rtol=0,
    )
    boundary_work = 0
    for face, cell, sign in ((1, 2, -1), (phi.sizes["x"] - 3, phi.sizes["x"] - 3, 1)):
        assert float(abs(faces["x"].isel(x=face)).max()) > 1e-6
        assert_allclose(
            first["x"].isel(x=face), weight.isel(x=cell) * faces["x"].isel(x=face)
        )
        boundary_work += sign * float(
            (
                faces["x"].isel(x=face) * (Phi.isel(x=face) + Phi.isel(x=face + 1)) / 2
            ).sum()
        )
    integral = lambda f: float((f * volume).isel(x=slice(2, -2)).sum())
    assert_allclose(
        integral(face_work(faces, Phi, g)),
        boundary_work - integral(Phi * face_divergence(faces, g)),
        rtol=1e-13,
        atol=1e-13,
    )
    complement = weighted_faces(faces, 1 - weight)
    assert_allclose(
        face_divergence({d: first[d] + complement[d] for d in faces}, g).isel(
            x=slice(2, -2)
        ),
        face_divergence(faces, g).isel(x=slice(2, -2)),
        atol=1e-13 * scale,
    )


def test_weighted_boundary_divergence_remains_bounded():
    errors = []
    for n in (16, 32, 64):
        v, g = sample("annulus", n)
        wall = SimpleNamespace(mxg=2)
        phi = fill_wall_guards(v["phi"], wall, "zero")
        psi = fill_wall_guards(v["psi"], wall, "pressure")
        omega = fill_wall_guards(v["omega"], wall, "neumann")
        faces = reynolds_face_fluxes(v["c"], omega, phi, psi, g)
        transport = face_divergence(weighted_faces(faces, 1 + 0.1 * np.sin(phi)), g)
        errors.append(float(abs(transport).isel(x=[2, -3]).mean()))
    # The full wall current remains bounded as dx decreases. This is a
    # boundedness check; the separate curvature test checks truncation orders.
    assert max(errors) < 2 * errors[0], errors


@pytest.mark.parametrize("n", (32, 64))
def test_upwind_omega_dissipates_periodic_kinetic_energy(n):
    v, g = sample("periodic", n)
    for name, value in {"J": 1.0, "g11": 1.0, "g33": 1.0, "g13": 0.0}.items():
        g[name] = value
    phi = periodic_halos(v["phi"])
    c = xr.ones_like(phi)
    omega = periodic_halos(perpendicular_operator(c, phi, g).divergence)
    faces = reynolds_face_fluxes(
        c,
        omega,
        phi,
        xr.zeros_like(phi),
        g,
        periodic_x=True,
    )
    work = (face_work(faces, phi, g) * g.J * g.dx * g.dy * g.dz).isel(x=slice(2, -2))
    assert float(work.sum()) > 0


def wall_curvature_state(
    n, nz=None, *, recovered=True, analytic_guards=False, zero_wall_curvature=False
):
    """FV-recovered flat-wall state with nonzero normal Reynolds current."""
    nz = n if nz is None else nz
    h = 1 / n
    x = (np.arange(n + 4) - 1.5) * h
    z = 2 * np.pi * np.arange(nz) / nz
    radial = x[:, None] * (1 - x[:, None])
    if zero_wall_curvature:
        # Preserve f'(0)=1, f'(1)=-1 and Omega_x=0, hence the same
        # nonzero wall current, while removing f'' at both physical walls.
        rx = 1 - 2 * x[:, None]
        fr = 1 + 2 * radial + 6.5 * radial**2
        frr = 2 + 13 * radial
        f = radial + radial**2 + (13 / 6) * radial**3
        fp = fr * rx
        fpp = frr * rx**2 - 2 * fr
        fppp = 13 * rx**3 - 6 * frr * rx
    else:
        f = radial - radial**2 / 12
        fp = (1 - 2 * x[:, None]) * (1 - radial / 6)
        fpp = -2 + radial / 3 - (1 - 2 * x[:, None]) ** 2 / 6
        fppp = 1 - 2 * x[:, None]
    # Recover Phi from the same FV Dirichlet problem as production. Sampling
    # the continuum potential alone would test an inconsistent boundary state.
    kz = 4 * np.sin(np.pi / nz) ** 2 / (2 * np.pi / nz) ** 2
    matrix = (
        np.diag(np.full(n, -2 / h**2 - kz))
        + np.diag(np.full(n - 1, 1 / h**2), 1)
        + np.diag(np.full(n - 1, 1 / h**2), -1)
    )
    matrix[0, 0] = matrix[-1, -1] = -3 / h**2 - kz
    generalised = np.zeros_like(f)
    generalised[2:-2, 0] = np.linalg.solve(matrix, (fpp - f)[2:-2, 0])
    if not recovered:
        generalised = f.copy()
    phi = xr.DataArray(
        generalised * (np.cos(z) + np.sin(z)) - f * np.sin(z), dims=("x", "z")
    )
    psi = xr.DataArray(1 + f * np.sin(z), dims=phi.dims)
    vort = xr.DataArray((fpp - f) * (np.cos(z) + np.sin(z)), dims=phi.dims)
    g = xr.Dataset(
        {
            name: xr.ones_like(phi) * value
            for name, value in {
                "J": 1,
                "dx": h,
                "dy": 1,
                "dz": 2 * np.pi / nz,
                "g11": 1,
                "g33": 1,
                "g13": 0,
            }.items()
        }
    )
    wall = SimpleNamespace(mxg=2)
    if not analytic_guards:
        phi = fill_wall_guards(phi, wall, "zero")
        psi = fill_wall_guards(psi, wall, "pressure")
        vort = fill_wall_guards(vort, wall, "neumann")
    # Direct differentiation of (v_E dot grad) grad(phi+psi).
    expected = xr.DataArray(
        (f * fppp - fp * fpp) / 2 * (np.cos(2 * z) - np.sin(2 * z))
        - (3 * fp * fpp + f * fppp) / 2,
        dims=phi.dims,
    )

    return phi, psi, vort, g, expected, radial, fpp, f, z


def wall_curvature_errors(n):
    phi, psi, vort, g, expected, radial, fpp, f, z = wall_curvature_state(n)
    h = 1 / n
    wall = SimpleNamespace(mxg=2)
    faces = reynolds_face_fluxes(1, vort, phi, psi, g)

    def norms(error):
        error = np.asarray(error)[2:-2]
        return {
            "boundary": float(np.max(abs(error[[0, -1]]))),
            "outer_three": float(np.sqrt(np.mean(error[[0, 1, 2, -3, -2, -1]] ** 2))),
            "interior": float(np.sqrt(np.mean(error[3:-3] ** 2))),
            "all": float(np.mean(abs(error))),
        }

    reynolds = face_divergence(faces, g)
    viscous = fv_xz(1, vort, g, boundary_flux=False).divergence
    viscous_exact = xr.DataArray(
        (-2 - 2 * fpp + f) * (np.cos(z) + np.sin(z)), dims=phi.dims
    )
    # A pure DC Dirichlet chi has nonzero wall current even with AC Neumann.
    # This is the exact FV Dirichlet solution of L(chi)=-2; sampling x(1-x)
    # alone would have a boundary truncation error from the midpoint reflection.
    chi = fill_wall_guards(
        xr.DataArray(radial[:, 0] + h * h / 4, dims="x") + xr.zeros_like(phi),
        wall,
        "zero",
    )
    chi_div = fv_xz(1, chi, g, boundary_flux=True).divergence
    return {
        "reynolds": norms(reynolds - expected),
        "viscous": norms(viscous - viscous_exact),
        "density_source": norms(reynolds + chi_div - (expected - 2)),
    }


def test_viscosity_with_wall_curvature_converges():
    errors = [wall_curvature_errors(n) for n in (32, 64, 128)]
    for region in ("boundary", "outer_three", "all", "interior"):
        values = np.asarray([e["viscous"][region] for e in errors])
        assert np.all(np.log2(values[:-1] / values[1:]) > 0.9), (
            region,
            errors,
        )


@pytest.mark.parametrize("term", ("reynolds", "density_source"))
def test_reynolds_curved_wall_limits_global_convergence(term):
    """O(1) first-cell error occupies O(h) volume; domain MAE is first order.

    Retain the original recovered state and nonzero wall current. The smooth
    interior converges near second order, but reflection of nonzero wall
    curvature leaves a bounded, nonvanishing boundary error. RMS over a fixed
    number of wall cells therefore also remains bounded rather than vanishing.
    """
    errors = [wall_curvature_errors(n)[term] for n in (32, 64, 128)]
    for region, lower, upper in (("all", 0.9, 1.5), ("interior", 1.7, 2.2)):
        values = np.asarray([e[region] for e in errors])
        orders = np.log2(values[:-1] / values[1:])
        assert np.all((orders > lower) & (orders < upper)), (region, errors, orders)
    for region in ("boundary", "outer_three"):
        values = np.asarray([e[region] for e in errors])
        assert np.all((values > 0.5 * values[0]) & (values < 2 * values[0])), (
            region,
            errors,
        )
    print("nonzero wall curvature", term, errors)


@pytest.mark.parametrize("limiter", ("mc", "central"))
def test_reynolds_zero_wall_curvature_converges(limiter):
    """Reflection-compatible wall accuracy, retaining FV recovery and MC.

    Wall maxima and all-cell mean absolute error should approach second order.
    MC-limited extrema retain lower RMS/maximum orders in the interior. A
    centred-donor diagnostic control checks second order in those norms too.
    """
    v = reference.reference("wall_zero_curvature")
    x, z = reference.x, reference.z
    for wall in (0, 1):
        assert s.simplify(v["phi"].subs(x, wall)) == 0
        assert s.simplify(v["psi"].subs(x, wall)) == 1
        for name in ("phi", "psi"):
            assert s.simplify(s.diff(v[name], x, 2).subs(x, wall)) == 0
        assert s.simplify(s.diff(v["omega"], x).subs(x, wall)) == 0
        assert s.simplify(v["expected_x"].subs(x, wall)) != 0

    exact = s.lambdify((x, z), v["expected"], "numpy", cse=True)
    errors = []
    for n in (32, 64, 128):
        phi, psi, omega, g, expected, *_ = wall_curvature_state(
            n, zero_wall_curvature=True
        )
        xx = (np.arange(n + 4) - 1.5) / n
        zz = 2 * np.pi * np.arange(n) / n
        # Check the closed-form divergence against independent covariant
        # differentiation, not a reference formed with numerical operators.
        assert_allclose(
            expected, exact(xx[:, None], zz[None, :]), rtol=1e-12, atol=1e-12
        )
        faces = reynolds_face_fluxes(1, omega, phi, psi, g, limiter=limiter)
        error = np.asarray(face_divergence(faces, g) - expected)[2:-2]
        errors.append(
            {
                "boundary": float(abs(error[[0, -1]]).max()),
                "outer_three": float(
                    np.sqrt(np.mean(error[[0, 1, 2, -3, -2, -1]] ** 2))
                ),
                "all": float(np.mean(abs(error))),
                "interior": float(np.sqrt(np.mean(error[3:-3] ** 2))),
                "rms": float(np.sqrt(np.mean(error**2))),
                "max": float(abs(error).max()),
            }
        )
        mismatch = omega - fv_xz(1, phi + psi, g).divergence
        assert float(abs(mismatch).isel(x=slice(2, -2)).max()) <= 1e-10
        # The wall current remains nontrivial after numerical recovery.
        assert float(abs(faces["x"].isel(x=[1, -3]) / g.dz.isel(x=[1, -3])).max()) > 0.1
        print("zero wall curvature", limiter, n, errors[-1])
    regions = (
        ("boundary", "outer_three", "all")
        if limiter == "mc"
        else ("boundary", "all", "rms", "max")
    )
    for region in regions:
        values = np.asarray([e[region] for e in errors])
        orders = np.log2(values[:-1] / values[1:])
        assert np.all(orders > 1.8), (region, errors, orders)


def test_flat_wall_covariant_reference_identities():
    """Verify tensor/sign/wall identities independently of compact R stencils."""
    v = reference.reference("wall")
    x, z = reference.x, reference.z
    Phi = v["phi"] + v["psi"]
    Rx = -s.diff(v["phi"], z) * s.diff(Phi, z, 2) - s.diff(v["phi"], x) * s.diff(
        Phi, x, z
    )
    Rz = s.diff(v["phi"], x) * s.diff(Phi, x, 2) + s.diff(v["phi"], z) * s.diff(
        Phi, x, z
    )
    assert (
        s.trigsimp(s.expand(v["expected_x"] - v["omega"] * s.diff(v["phi"], z) - Rx))
        == 0
    )
    assert (
        s.trigsimp(s.expand(v["expected_z"] + v["omega"] * s.diff(v["phi"], x) - Rz))
        == 0
    )
    for wall in (0, 1):
        assert s.simplify(s.diff(v["omega"], x).subs(x, wall)) == 0
        assert s.simplify(v["expected_x"].subs(x, wall)) != 0


@pytest.mark.parametrize("recovered", (False, True))
def test_reflected_wall_leading_error(recovered):
    """Verify the accepted bounded wall error against its independent limit.

    Independent Taylor expansion: f''(0)=-13/6. FV-recovered Phi has a
    -f''*h^2/8 boundary shift, while sampled psi retains f''*h^2/8. Separate
    reflection therefore creates an O(h) wall face discrepancy. Dividing by
    the first cell width leaves the nonzero limits checked here.
    """
    differences = []
    for n in (32, 64, 128):
        phi, psi, omega, g, expected, *_ = wall_curvature_state(n, recovered=recovered)
        divergence = face_divergence(reynolds_face_fluxes(1, omega, phi, psi, g), g)
        z = 2 * np.pi * np.arange(n) / n
        leading = (
            13 / 24 * np.sin(z) * (np.cos(z) - np.sin(z))
            if recovered
            else np.full(n, -13 / 24)
        )
        error = (divergence - expected).isel(x=2).values
        differences.append(float(np.max(abs(error - leading))))
    # The residual approaches the Taylor limit; this does not reclassify the
    # finite, nonzero wall error as a passing convergence result.
    assert np.all(np.log2(np.asarray(differences[:-1]) / differences[1:]) > 0.9), (
        differences
    )
