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


def chi_shape(a, theta):
    """Analytic cosine particular satisfying both X face conditions.

    For phi=sin(a)*cos(theta), psi=P(a)*sin(theta), kZ/kX=2:
    (d_aa+4*d_theta_theta)chi = physical Reynolds divergence / c.
    Multiply by kX^2 and the two amplitudes for physical chi.
    """
    return (
        40 * s.cos(5 * a)
        - s.Rational(792, 7) * s.cos(3 * a)
        + s.Rational(512, 7) * s.cos(a)
    ) / 41 + (84 * s.cos(5 * a) - 252 * s.cos(3 * a) + 168 * s.cos(a)) * s.sin(
        theta
    ) ** 2 / 41


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
    ni = 1 + r(5, 100) * s.cos(2 * a) * s.cos(b - tau) + r(3, 100) * s.cos(d + tau)
    nj = r(1, 4) * (1 + r(5, 100) * s.cos(2 * a) * s.sin(theta) * s.cos(b + tau))
    vi = r(1, 10) + r(2, 100) * s.cos(2 * a) * s.sin(b - tau) * s.cos(theta)
    vj = r(7, 100) + r(15, 1000) * s.cos(2 * a) * s.cos(b + tau) * s.sin(theta)
    te = 1 + r(5, 100) * s.cos(2 * a) * s.sin(theta) * s.cos(b - tau)
    amp_phi = r(2, 10) * (1 + r(1, 10) * s.sin(b + tau))
    amp_psi = r(2, 100) * (1 + r(5, 100) * s.cos(b + tau))
    phi = amp_phi * s.sin(a) * s.cos(theta)
    shape = r(813, 287) - r(1100, 287) * s.cos(2 * a) + s.cos(4 * a)
    psi = r(12, 10) + r(4, 100) * s.cos(b - tau) + amp_psi * shape * s.sin(theta)
    if variant != "baseline":
        # These focused recovery/viscosity cases vanish to sixth order at X
        # faces: all mixed normal fluxes and Omega_X are compatible there.
        phi = s.Integer(0)
        psi = (
            r(12, 10)
            + r(4, 100) * s.cos(b - tau)
            + amp_psi * (1 - s.cos(2 * a)) ** 3 * s.sin(theta)
        )
    ti = common_temperature(psi, ni, nj)
    pi, pj = ni * ti, nj * ti
    ne = ni + 2 * nj
    ji = ni * vi + 2 * nj * vj
    ve, je = ji / ne, -ji
    assert s.cancel(ji - ne * ve) == 0, "Zero total parallel current"
    pe = ne * te
    mi, mj = 2 * ni * vi, 3 * nj * vj
    ki, kj = ni * vi**2, r(3, 2) * nj * vj**2
    rho, pressure, Phi, c = 2 * ni + 3 * nj, pi + pj, phi + psi, s.Integer(2)
    lap_phi = -5 * (s.pi / LENGTH) ** 2 * phi
    omega = c * (lap_phi + lap(psi, H))
    kperp = c * grad(Phi).dot(H * grad(Phi)) / 2
    energy = r(3, 2) * pressure + ki + kj + kperp
    alpha = advection(psi, phi)
    raw_reynolds = (
        -omega * drift(phi) / 2 - c * grad(alpha) / 2 - c * lap_phi * drift(Phi) / 2
    )
    # Combine the analytic Reynolds terms before exporting them. This removes
    # cancelling products without changing their continuum or boundary limit.
    angular = s.Symbol("angular", real=True)
    poly = r(813, 287) - r(1100, 287) * s.cos(2 * angular) + s.cos(4 * angular)
    prime = s.diff(poly, angular).subs(angular, a)
    second = s.diff(poly, angular, 2).subs(angular, a)
    k = s.pi / LENGTH
    reynolds = (
        c
        * k**3
        * s.Matrix(
            [
                -10 * amp_phi**2 * s.sin(a) ** 2 * s.cos(theta) * s.sin(theta)
                + amp_phi
                * amp_psi
                * (
                    prime * s.cos(a)
                    + (2 * second - 4 * shape) * s.sin(a) * s.sin(theta) ** 2
                    + 4 * shape * s.sin(a) * s.cos(theta) ** 2
                ),
                -5 * amp_phi**2 * s.sin(a) * s.cos(a) * s.cos(theta) ** 2
                + amp_phi
                * amp_psi
                * ((second / 2 - 6 * shape) * s.cos(a) + r(3, 2) * prime * s.sin(a))
                * s.sin(theta)
                * s.cos(theta),
            ]
        )
    )
    test_point = {X: 1.37, Y: 2.91, Z: 4.23, t: 0.73}
    if variant == "baseline":
        assert (
            max(abs(float(q.subs(test_point))) for q in reynolds - raw_reynolds) < 1e-14
        )
    chi = (s.pi / LENGTH) ** 2 * amp_phi * amp_psi * chi_shape(a, theta)
    if variant != "baseline":
        reynolds = s.zeros(2, 1)
        chi = nu * omega / c
    jnu = nu * H * grad(omega)
    jpol = -c * H * grad(chi) + reynolds
    terms = {"Vort.reynolds": div(reynolds)}
    # The verified finite cosine inversion makes these expressions identical.
    # Keeping the harmonic form avoids large cancelling terms in the source.
    rhs = {"Vort": c * lap(chi, H)}
    if nu:
        terms["Vort.viscosity"] = div(jnu)

    def polarisation(weight):
        return -grad(weight).dot(jpol) if not nu else -div(weight * jpol)

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
            # Analytic div(Jpol)=0, as verified by the auxiliary-current identity.
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
    a, d = s.symbols("a d", real=True)
    shape = s.Rational(813, 287) - s.Rational(1100, 287) * s.cos(2 * a) + s.cos(4 * a)
    f, p = s.sin(a) * s.cos(d), shape * s.sin(d)
    L = lambda q: s.diff(q, a, 2) + 4 * s.diff(q, d, 2)
    D = lambda q, h: 2 * (s.diff(q, a) * s.diff(h, d) - s.diff(q, d) * s.diff(h, a))
    alpha = D(p, f)
    forcing = -(D(L(f + p), f) + L(alpha) + D(L(f), f + p)) / 2
    chi = chi_shape(a, d)
    residual = s.expand(s.expand_trig(L(chi) - forcing).rewrite(s.exp))
    assert s.simplify(residual) == 0, "Analytic chi does not solve the physical RHS"
    rx = -L(f + p) * s.diff(f, d) - s.diff(alpha, a) / 2 - L(f) * s.diff(f + p, d)
    for face in (0, s.pi):
        for q in (shape, s.diff(shape, a), chi, s.diff(chi, a), rx, s.diff(alpha, a)):
            assert s.simplify(q.subs(a, face)) == 0
    fields, rhs, _, terms = model(variant)
    H = metric(variant)
    if variant != "baseline":
        for face in (0, LENGTH):
            for name in ("psi", "Vort", "chi"):
                assert s.simplify(s.diff(fields[name], X).subs(X, face)) == 0
            for name in ("chi", "Vort"):
                assert s.simplify(fields[name].subs(X, face)) == 0
            assert s.simplify((H * grad(fields["psi"]))[0].subs(X, face)) == 0
        assert s.simplify(2 * lap(fields["chi"], H) - rhs["Vort"]) == 0
    for face in (0, LENGTH):
        assert s.simplify(fields["phi"].subs(X, face)) == 0
        for name in ("Ni", "Nj", "NVi", "NVj", "Pe", "Pi", "pi", "pj", "Tion", "psi"):
            assert s.simplify(s.diff(fields[name], X).subs(X, face)) == 0, name
    # Rigorous primitive bounds (|P| <= 2200/287); density floors stay inactive.
    assert (
        s.Rational(12, 10)
        - s.Rational(4, 100)
        - s.Rational(2, 100) * s.Rational(2200, 287) * s.Rational(105, 100)
        > 0
    )
    sample = (1.37, 2.91, 4.23, 0.73)
    magnitudes = {
        k: abs(float(s.lambdify(COORDINATES, v, "numpy", cse=True)(*sample)))
        for k, v in terms.items()
    }
    assert all(np.isfinite(v) for v in magnitudes.values()), magnitudes
    if variant == "baseline":
        assert all(v > 1e-12 for v in magnitudes.values()), magnitudes
    else:
        for key in ("Ni.parallel", "NVj.pressure", "Pe.compression", "Eion.parallel"):
            assert magnitudes[key] > 1e-12, key
        if variant == "viscosity":
            for key in (
                "Vort.viscosity",
                "Ni.polarisation",
                "Eion.viscous_flux",
                "Pe.viscous_heating",
            ):
                assert magnitudes[key] > 1e-12, key
    return {
        "variant": variant,
        "elliptic_and_X_boundaries": "verified symbolically",
        "positive_primitives": True,
        "sample_term_magnitudes": magnitudes,
    }
