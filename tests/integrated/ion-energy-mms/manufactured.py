"""Continuum two-ion slab MMS, differentiated in normalised physical coordinates.

No numerical operator or finest-grid solution enters the reference or sources.
The auxiliary chi solves the physical current equation before BOUT MMS forcing.
"""

from functools import cache, lru_cache

import numpy as np
import sympy as s

X, Y, Z, t = s.symbols("X Y Z t", real=True)
COORDINATES = (X, Y, Z, t)
LENGTH = s.Integer(10)
EVOLVED = ("Ni", "Nj", "NVi", "NVj", "Pe", "Vort", "Eion")
RECOVERED = ("Pi", "pi", "pj", "Tion", "phi", "psi", "chi")
VARIANTS = ("baseline", "skew", "viscosity")


def grad(f):
    return s.Matrix([s.diff(f, X), s.diff(f, Z)])


def div(v):
    return s.diff(v[0], X) + s.diff(v[1], Z)


def metric(variant="baseline"):
    sigma = s.Rational(1, 4) if variant == "skew" else s.Integer(0)
    return s.Matrix([[1, sigma], [sigma, 1 + sigma**2]])


def lap(f, H=None):
    return div((s.eye(2) if H is None else H) * grad(f))


def drift(f):
    return s.Matrix([s.diff(f, Z), -s.diff(f, X)])


def advection(q, f):
    return grad(q).dot(drift(f))


def common_temperature(psi, ni, nj):
    return 2 * psi / (2 * ni + s.Rational(3, 2) * nj)


@lru_cache(maxsize=3)
def model(variant="baseline"):
    """Return exact fields, physical RHS, forcing and separate term coverage."""
    r = s.Rational
    a, b, d, tau = (
        s.pi * X / LENGTH,
        2 * s.pi * Y / LENGTH,
        2 * s.pi * Z / LENGTH,
        t / 5,
    )
    theta = d - tau
    if variant not in VARIANTS:
        raise ValueError(f"Unknown MMS variant: {variant}")
    H = metric(variant)
    nu = r(2, 100) if variant == "viscosity" else s.Integer(0)
    shape = (1 - s.cos(2 * a)) ** 2
    ni = 1 + r(5, 100) * s.cos(2 * a) + r(3, 400) * shape * s.cos(theta)
    nj = ni / 4
    # Small zero-mean parallel flows keep acoustic gradients active without
    # making limiter errors at advected extrema dominate continuum convergence.
    vi = r(1, 1000) * (s.sin(b - tau) + s.sin(3 * (b - tau)) / 9)
    vj = r(1, 1000) * (s.cos(b + tau) - s.cos(3 * (b + tau)) / 9)
    phi = s.Integer(0)
    psi_mean = r(12, 10) + r(4, 1000) * s.cos(b - tau)
    psi_fluctuation = r(2, 100) * shape * s.sin(theta)
    psi = psi_mean + psi_fluctuation
    ti = common_temperature(psi, ni, nj)
    pi, pj = ni * ti, nj * ti
    ne = ni + 2 * nj
    ji = ni * vi + 2 * nj * vj
    ve, je = ji / ne, -ji
    # Prescribe a closed curvature-driven auxiliary flow, then obtain the
    # electron pressure from c Lap(chi)=curl_z d_Z(Pe+Pi). The antiderivative
    # is analytic and both chi and its normal derivative vanish at each wall.
    curl_z = r(1, 10)
    chi_curvature = r(1, 100) * shape * s.cos(theta)
    wave = 2 * s.pi / LENGTH
    total_pressure_fluctuation = (
        r(2, 1)
        / curl_z
        * r(1, 100)
        * (
            (s.diff(shape, X, 2) - H[1, 1] * wave**2 * shape) * s.sin(theta) / wave
            + 2 * H[0, 1] * s.diff(shape, X) * s.cos(theta)
        )
    )
    pe = 4 - r(20, 19) * psi_fluctuation + total_pressure_fluctuation
    mi, mj = 2 * ni * vi, 3 * nj * vj
    ki, kj = ni * vi**2, r(3, 2) * nj * vj**2
    rho, pressure, Phi, c = 2 * ni + 3 * nj, pi + pj, phi + psi, s.Integer(2)
    lap_phi = lap(phi, H)
    omega = c * (lap_phi + lap(psi, H))
    kperp = c * grad(Phi).dot(H * grad(Phi)) / 2
    energy = r(3, 2) * pressure + ki + kj + kperp
    alpha = s.Integer(0)
    reynolds = s.zeros(2, 1)
    chi = chi_curvature + nu * omega / c
    jnu = nu * H * grad(omega)
    jpol = -c * H * grad(chi) + reynolds
    terms = {"Vort.curvature": curl_z * s.diff(pe + pressure, Z)}
    # Use the verified auxiliary identity to avoid cancelling source terms.
    rhs = {"Vort": c * lap(chi, H)}
    if nu:
        terms["Vort.viscosity"] = div(jnu)

    def polarisation(weight):
        return -div(weight * jpol)

    for name, n, v, m, p, mass, charge in (
        ("i", ni, vi, mi, pi, 2, 1),
        ("j", nj, vj, mj, pj, 3, 2),
    ):
        w = mass / (charge * rho)
        terms.update(
            {
                f"N{name}.exb": -advection(n, phi),
                f"N{name}.parallel": -s.diff(n * v, Y),
                f"N{name}.polarisation": polarisation(w * n),
                f"NV{name}.exb": -advection(m, phi),
                f"NV{name}.parallel": -s.diff(m * v, Y),
                f"NV{name}.pressure": -s.diff(p, Y),
                f"NV{name}.electric": -charge * n * s.diff(phi, Y),
                f"NV{name}.polarisation": polarisation(w * m),
            }
        )
        for variable in (f"N{name}", f"NV{name}"):
            rhs[variable] = sum(
                v for k, v in terms.items() if k.startswith(variable + ".")
            )
    terms.update(
        {
            "Pe.exb": -advection(pe, phi),
            "Pe.parallel": -s.diff(pe * ve, Y),
            "Pe.compression": -r(2, 3) * pe * s.diff(ve, Y),
            "Eion.exb": -advection(r(3, 2) * pressure + ki + kj, phi),
            "Eion.reynolds": -div(c * alpha * grad(Phi) / 2),
            "Eion.parallel_electric": -s.diff(phi * ji, Y),
            "Eion.parallel": -s.diff(
                (r(5, 2) * pi + ki) * vi + (r(5, 2) * pj + kj) * vj, Y
            ),
            # The auxiliary identity gives div(Jpol) = -physical ddt(Vort).
            "Eion.polarisation": polarisation(
                phi
                + 2 / rho * (r(5, 2) * pi + ki)
                + r(3, 2) / rho * (r(5, 2) * pj + kj)
            ),
            "Eion.electron_work": -phi * s.diff(je, Y),
        }
    )
    if nu:
        heating = -nu * grad(omega).dot(H * grad(Phi))
        electron_heating = (ne / 1836) / (rho + ne / 1836) * heating
        terms.update(
            {
                "Eion.viscous_flux": -div(Phi * jnu),
                "Eion.electron_heating": -electron_heating,
                "Pe.viscous_heating": r(2, 3) * electron_heating,
            }
        )
    for variable in ("Pe", "Eion"):
        rhs[variable] = sum(v for k, v in terms.items() if k.startswith(variable + "."))
    fields = dict(zip(EVOLVED, (ni, nj, mi, mj, pe, omega, energy)))
    fields.update(Pi=pressure, pi=pi, pj=pj, Tion=ti, phi=phi, psi=psi, chi=chi)
    sources = {name: s.diff(fields[name], t) - rhs[name] for name in EVOLVED}
    return fields, rhs, sources, terms


@cache
def evaluator(name, variant="baseline"):
    return s.lambdify(COORDINATES, model(variant)[0][name], "numpy", cse=True)


def evaluate(name, x, y, z, time, variant="baseline"):
    """Evaluate on broadcastable physical coordinates, including guards."""
    result = evaluator(name, variant)(x, y, z, time)
    return np.broadcast_to(
        result,
        np.broadcast_shapes(np.shape(x), np.shape(y), np.shape(z), np.shape(time)),
    )


def preflight(variant="baseline"):
    """Check the analytic elliptic, closure and boundary contracts."""
    fields, rhs, _, terms = model(variant)
    H = metric(variant)
    assert s.trigsimp(2 * lap(fields["chi"], H) - rhs["Vort"]) == 0
    residual = 2 * lap(fields["chi"], H) - sum(
        value for key, value in terms.items() if key.startswith("Vort.")
    )
    assert s.simplify(s.expand(residual)) == 0, "Auxiliary current identity"
    for face in (0, LENGTH):
        for name in (*EVOLVED, *RECOVERED):
            assert s.simplify(s.diff(fields[name], X).subs(X, face)) == 0, name
        for name in ("phi", "Vort", "chi"):
            assert s.simplify(fields[name].subs(X, face)) == 0, name
        assert s.simplify((H * grad(fields["chi"]))[0].subs(X, face)) == 0
    # Analytic amplitude bounds for density, psi and electron pressure.
    psi_amplitude = s.Rational(2, 100)
    assert 1 - s.Rational(5, 100) - s.Rational(3, 400) * 4 > 0
    assert s.Rational(12, 10) - s.Rational(4, 1000) - psi_amplitude * 4 > 0
    assert (
        4 - s.Rational(80, 19) * psi_amplitude - s.Rational(39, 10) * s.pi / LENGTH > 0
    )
    sample = (1.37, 2.91, 4.23, 0.73)
    magnitudes = {
        k: abs(float(s.lambdify(COORDINATES, v, "numpy", cse=True)(*sample)))
        for k, v in terms.items()
    }
    assert all(np.isfinite(v) for v in magnitudes.values()), magnitudes
    for key in (
        "Vort.curvature",
        "Ni.parallel",
        "Ni.polarisation",
        "NVj.pressure",
        "Pe.compression",
        "Eion.parallel",
        "Eion.polarisation",
    ):
        assert magnitudes[key] > 1e-12, key
    if variant == "viscosity":
        for key in ("Vort.viscosity", "Eion.viscous_flux", "Pe.viscous_heating"):
            assert magnitudes[key] > 1e-12, key
    return {
        "variant": variant,
        "elliptic_and_X_boundaries": "verified symbolically",
        "positive_primitives": True,
        "sample_term_magnitudes": magnitudes,
    }
