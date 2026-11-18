"""Continuum two-ion slab MMS, differentiated in normalised physical coordinates.

No numerical operator or finest-grid solution enters the reference or sources.
The auxiliary chi excludes viscous stress, which has no particle flux,
before BOUT MMS forcing. Viscosity still diffuses vorticity and exchanges energy.
"""

from functools import cache, lru_cache

import numpy as np
import sympy as s
from cases import EVOLVED, PARALLEL_VARIANTS, RECOVERED, VARIANTS

X, Y, Z, t = s.symbols("X Y Z t", real=True)
COORDINATES = (X, Y, Z, t)
LENGTH = s.Integer(10)


def evolved(variant):
    return EVOLVED + (("NVe",) if variant in PARALLEL_VARIANTS else ())


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


def pressure_shape(a, variant):
    if variant == "parallel_current":
        return 1 + (1 - s.cos(2 * a)) ** 2
    if variant == "parallel_current_wall":
        return 2 - s.cos(2 * a)
    return (
        1 + s.Rational(813, 287) - s.Rational(1100, 287) * s.cos(2 * a) + s.cos(4 * a)
    )


def chi_shape(a, theta, *, variant="baseline"):
    """Analytic cosine particular satisfying both X face conditions.

    For phi=sin(a)*cos(theta), psi=P(a)*sin(theta), kZ/kX=2:
    (d_aa+4*d_theta_theta)chi = physical Reynolds divergence / c.
    Multiply by kX^2 and the two amplitudes for physical chi.
    """
    wall_harmonic = s.Rational(8, 17) * s.cos(a) * s.cos(2 * theta)
    if variant == "parallel_current":
        # The constant radial pressure harmonic has nonzero wall AC. The
        # remaining shape has P'=P''=0 at both radial walls.
        # The cosine inverse has Neumann AC and zero Dirichlet DC traces.
        return (
            wall_harmonic
            + s.Rational(16, 17) * s.cos(a)
            - s.Rational(36, 25) * s.cos(3 * a)
            + s.Rational(20, 41) * s.cos(5 * a)
            + (
                s.Rational(36, 17) * s.cos(a)
                - s.Rational(78, 25) * s.cos(3 * a)
                + s.Rational(42, 41) * s.cos(5 * a)
            )
            * s.sin(theta) ** 2
        )
    if variant == "parallel_current_wall":
        # P=2-cos(2a) retains finite pressure-potential curvature at both
        # walls, with the same nonzero wall AC as the zero-curvature profile.
        return (
            wall_harmonic
            + s.Rational(10, 17) * s.cos(a)
            - s.Rational(18, 25) * s.cos(3 * a)
            + (s.Rational(14, 17) * s.cos(a) - s.Rational(14, 25) * s.cos(3 * a))
            * s.sin(theta) ** 2
        )
    return (
        wall_harmonic
        + (
            40 * s.cos(5 * a)
            - s.Rational(792, 7) * s.cos(3 * a)
            + s.Rational(512, 7) * s.cos(a)
        )
        / 41
        + (84 * s.cos(5 * a) - 252 * s.cos(3 * a) + 168 * s.cos(a))
        * s.sin(theta) ** 2
        / 41
    )


@lru_cache(maxsize=len(VARIANTS))
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
    # Wall-uniform density and velocity isolate the parallel-current operator
    # comparison; temperature retains the local pressure-potential fluctuation.
    wall_envelope = s.sin(a) ** 2 if variant in PARALLEL_VARIANTS else s.Integer(1)
    ni = (
        1
        + r(5, 100) * s.cos(2 * a) * s.cos(b - tau)
        + r(3, 100) * wall_envelope * s.cos(d + tau)
    )
    nj = r(1, 4) * (
        1 + r(5, 100) * wall_envelope * s.cos(2 * a) * s.sin(theta) * s.cos(b + tau)
    )
    vi = r(1, 10) + r(2, 100) * wall_envelope * s.cos(2 * a) * s.sin(b - tau) * s.cos(
        theta
    )
    vj = r(7, 100) + r(15, 1000) * wall_envelope * s.cos(2 * a) * s.cos(
        b + tau
    ) * s.sin(theta)
    te = 1 + r(5, 100) * wall_envelope * s.cos(2 * a) * s.sin(theta) * s.cos(b - tau)
    amp_phi = r(2, 10) * (1 + r(1, 10) * s.sin(b + tau))
    amp_psi = r(2, 100) * (1 + r(5, 100) * s.cos(b + tau))
    phi = amp_phi * s.sin(a) * s.cos(theta)
    shape = pressure_shape(a, variant)
    psi = r(12, 10) + r(4, 100) * s.cos(b - tau) + amp_psi * shape * s.sin(theta)
    if variant in ("skew", "viscosity"):
        # The radial perturbation vanishes to sixth order at X faces; the
        # constant radial harmonic retains nonzero wall AC and Neumann guards.
        phi = s.Integer(0)
        psi = (
            r(12, 10)
            + r(4, 100) * s.cos(b - tau)
            + amp_psi * (1 + (1 - s.cos(2 * a)) ** 3) * s.sin(theta)
        )
    ti = common_temperature(psi, ni, nj)
    pi, pj = ni * ti, nj * ti
    ne = ni + 2 * nj
    ji = ni * vi + 2 * nj * vj
    ve, je = ji / ne, -ji
    assert s.cancel(ji - ne * ve) == 0, "Zero total parallel current"
    parallel_current = s.Integer(0)
    if variant in PARALLEL_VARIANTS:
        # An independent electron momentum leaves a finite current/divergence.
        # This Neumann X cosine and periodic Y/Z harmonic has a closed inverse.
        parallel_current = r(3, 100) * s.cos(2 * a) * s.sin(b - tau) * s.cos(theta)
        je = parallel_current - ji
        ve = -je / ne
    pe = ne * te
    curvature = r(1, 10) if variant == "curvature" else s.Integer(0)
    if curvature:
        # Cancel the composition-dependent pressure only in this added variant.
        # Total pressure is then a finite Neumann cosine mode with an analytic
        # inverse for the diamagnetic current; electron pressure stays positive.
        pe = 4 + r(1, 50) * s.cos(2 * a) * s.sin(theta) - pi - pj
    mi, mj = 2 * ni * vi, 3 * nj * vj
    ki, kj = ni * vi**2, r(3, 2) * nj * vj**2
    rho, pressure, Phi, c = 2 * ni + 3 * nj, pi + pj, phi + psi, s.Integer(2)
    lap_phi = -5 * (s.pi / LENGTH) ** 2 * phi
    omega = c * (lap_phi + lap(psi, H))
    kperp = c * grad(Phi).dot(H * grad(Phi)) / 2
    energy = r(3, 2) * pressure + ki + kj + kperp
    # Combine the analytic current harmonics before exporting them. This removes
    # cancelling products without changing their continuum or boundary limit.
    angular = s.Symbol("angular", real=True)
    poly = pressure_shape(angular, variant)
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
    chi = (
        (s.pi / LENGTH) ** 2 * amp_phi * amp_psi * chi_shape(a, theta, variant=variant)
    )
    if variant in ("skew", "viscosity"):
        reynolds = s.zeros(2, 1)
        chi = s.Integer(0)
    if curvature:
        chi += -curvature * r(1, 50) / (4 * c * k) * s.cos(2 * a) * s.cos(theta)
    if variant in PARALLEL_VARIANTS:
        # c*lap(chi_parallel) = d_Y(Jparallel), with kX=kZ=2*pi/L.
        chi += -r(3, 100) / (4 * c * k) * s.cos(2 * a) * s.cos(b - tau) * s.cos(theta)
    # Retain the divergence-free contribution to the covariant drift current.
    # It contributes to particle, momentum and energy transport through walls.
    stream = -c * grad(phi).dot(H * grad(Phi)) / 2
    reynolds -= drift(stream)
    jnu = nu * H * grad(omega)
    jpol = -c * H * grad(chi) + reynolds
    terms = {"Vort.reynolds": div(reynolds)}
    if variant in PARALLEL_VARIANTS:
        terms["Vort.parallel_current"] = s.diff(parallel_current, Y)
    # The verified finite cosine inversion makes these expressions identical.
    # Keeping the harmonic form avoids large cancelling terms in the source.
    rhs = {"Vort": c * lap(chi, H)}
    if curvature:
        terms["Vort.curvature"] = curvature * s.diff(pe + pressure, Z)
    if nu:
        terms["Vort.viscosity"] = div(jnu)
        rhs["Vort"] += terms["Vort.viscosity"]

    def polarisation(weight):
        return (
            -grad(weight).dot(jpol)
            if not curvature and variant not in PARALLEL_VARIANTS
            else -div(weight * jpol)
        )

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
            "Eion.reynolds": -advection(kperp, phi),
            "Eion.parallel_electric": -s.diff(phi * ji, Y),
            "Eion.parallel": -s.diff(
                (r(5, 2) * pi + ki) * vi + (r(5, 2) * pj + kj) * vj, Y
            ),
            # Analytic div(Jpol)=0, as verified by the auxiliary-current identity.
            "Eion.polarisation": polarisation(
                Phi
                + 2 / rho * (r(3, 2) * pi + ki)
                + r(3, 2) / rho * (r(3, 2) * pj + kj)
            ),
            "Eion.electron_work": -phi * (s.diff(je, Y) + curvature * s.diff(pe, Z)),
        }
    )
    if curvature:
        terms["Pe.curvature_work"] = -r(2, 3) * curvature * pe * s.diff(phi, Z)
        terms["Eion.curvature_electric"] = -curvature * s.diff(phi * pressure, Z)
    if nu:
        heating = -nu * grad(omega).dot(H * grad(Phi))
        electron_heating = (ne / 1836) / (rho + ne / 1836) * heating
        terms.update(
            {
                "Eion.electron_heating": -electron_heating,
                "Pe.viscous_heating": r(2, 3) * electron_heating,
            }
        )
    for variable in ("Pe", "Eion"):
        rhs[variable] = sum(v for k, v in terms.items() if k.startswith(variable + "."))
    fields = dict(zip(EVOLVED, (ni, nj, mi, mj, pe, omega, energy)))
    fields.update(Pi=pressure, pi=pi, pj=pj, Tion=ti, phi=phi, psi=psi, chi=chi)
    if variant in PARALLEL_VARIANTS:
        me = -je / 1836
        terms.update(
            {
                "NVe.exb": -advection(me, phi),
                "NVe.parallel": -s.diff(me * ve, Y),
                "NVe.pressure": -s.diff(pe, Y),
                "NVe.electric": ne * s.diff(phi, Y),
            }
        )
        rhs["NVe"] = sum(v for key, v in terms.items() if key.startswith("NVe."))
        fields.update(
            NVe=me,
            Te=pe / ne,
            Jpar=parallel_current,
            DivJpar=s.diff(parallel_current, Y),
        )
    sources = {name: s.diff(fields[name], t) - rhs[name] for name in evolved(variant)}
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
    shape = pressure_shape(a, "baseline")
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
        assert s.simplify(shape.subs(a, face)) == 1
        assert s.simplify(s.integrate(chi.subs(a, face), (d, 0, 2 * s.pi))) == 0
        for q in (s.diff(shape, a), s.diff(chi, a), rx, s.diff(alpha, a)):
            assert s.simplify(q.subs(a, face)) == 0
    fields, rhs, _, terms = model(variant)
    assert all(name in fields for name in RECOVERED), "Missing recovered MMS field"
    H = metric(variant)
    Phi = fields["phi"] + fields["psi"]
    v = drift(fields["phi"])
    u = H * grad(Phi)
    skew = 2 * s.Matrix([grad(u[i]).dot(v) for i in range(2)])
    residual = div(-skew) - terms["Vort.reynolds"]
    assert s.simplify(s.expand(s.expand_trig(residual).rewrite(s.exp))) == 0, (
        "Reynolds skew identity"
    )
    assert (
        s.simplify(
            terms["Eion.reynolds"]
            + advection(grad(Phi).dot(H * grad(Phi)), fields["phi"])
        )
        == 0
    )
    if variant == "curvature":
        residual = 2 * lap(fields["chi"], H) - sum(
            terms[key] for key in ("Vort.reynolds", "Vort.curvature")
        )
        assert s.simplify(s.expand(s.expand_trig(residual).rewrite(s.exp))) == 0, (
            "Curvature chi inversion"
        )
    if variant in PARALLEL_VARIANTS:
        residual = 2 * lap(fields["chi"], H) - sum(
            terms[key] for key in ("Vort.reynolds", "Vort.parallel_current")
        )
        assert s.simplify(s.expand(s.expand_trig(residual).rewrite(s.exp))) == 0, (
            "Parallel-current chi inversion"
        )
        current = fields["NVi"] / 2 + 2 * fields["NVj"] / 3 - 1836 * fields["NVe"]
        assert s.simplify(current - fields["Jpar"]) == 0
        assert s.simplify(s.diff(current, Y) - fields["DivJpar"]) == 0
        for face in (0, LENGTH):
            for name in ("Ni", "Nj", "NVi", "NVj", "Pe"):
                assert s.simplify(s.diff(fields[name], Z).subs(X, face)) == 0, (
                    name,
                    face,
                )
            assert s.simplify(s.diff(fields["chi"], X).subs(X, face)) == 0
            assert s.simplify(s.diff(fields["NVe"], X).subs(X, face)) == 0
            # The added chi is an AC cosine; its Z average vanishes at each wall.
            assert (
                s.simplify(s.integrate(fields["chi"].subs(X, face), (Z, 0, LENGTH)))
                == 0
            )
        for face in (0, LENGTH):
            for name, potential in (
                ("phi", fields["phi"]),
                ("psi", fields["psi"]),
                ("Phi", Phi),
            ):
                wall_curvature = s.simplify(s.diff(potential, X, 2).subs(X, face))
                if variant == "parallel_current" or name == "phi":
                    assert wall_curvature == 0, (name, face)
                else:
                    assert wall_curvature != 0, (name, face)

    if variant in ("skew", "viscosity"):
        for face in (0, LENGTH):
            for name in ("psi", "Vort", "chi"):
                assert s.simplify(s.diff(fields[name], X).subs(X, face)) == 0
        assert fields["chi"] == 0
        assert s.simplify(rhs["Vort"] - terms.get("Vort.viscosity", 0)) == 0
    for face in (0, LENGTH):
        assert s.simplify(fields["phi"].subs(X, face)) == 0
        assert s.simplify(s.diff(fields["psi"], Z).subs(X, face)) != 0
        for name in ("Ni", "Nj", "NVi", "NVj", "Pe", "Pi", "pi", "pj", "Tion", "psi"):
            assert s.simplify(s.diff(fields[name], X).subs(X, face)) == 0, name
    # Rigorous primitive bounds (|P| <= 9); floors stay inactive.
    assert (
        s.Rational(12, 10)
        - s.Rational(4, 100)
        - s.Rational(2, 100) * s.Integer(9) * s.Rational(105, 100)
        > 0
    )
    if variant == "curvature":
        psi_upper = (
            s.Rational(12, 10)
            + s.Rational(4, 100)
            + s.Rational(2, 100) * s.Integer(9) * s.Rational(105, 100)
        )
        # Pi <= 4*psi/3 for positive ni,nj, so the replaced Pe is positive too.
        assert 4 - s.Rational(1, 50) - s.Rational(4, 3) * psi_upper > 0
    sample = (1.37, 2.91, 4.23, 0.73)
    magnitudes = {
        k: abs(float(s.lambdify(COORDINATES, v, "numpy", cse=True)(*sample)))
        for k, v in terms.items()
    }
    assert all(np.isfinite(v) for v in magnitudes.values()), magnitudes
    if variant in ("baseline", "curvature", *PARALLEL_VARIANTS):
        assert all(v > 1e-12 for v in magnitudes.values()), magnitudes
    else:
        for key in ("Ni.parallel", "NVj.pressure", "Pe.compression", "Eion.parallel"):
            assert magnitudes[key] > 1e-12, key
        if variant == "viscosity":
            for key in (
                "Vort.viscosity",
                "Eion.electron_heating",
                "Pe.viscous_heating",
            ):
                assert magnitudes[key] > 1e-12, key
            for key in (
                "Ni.polarisation",
                "Nj.polarisation",
                "NVi.polarisation",
                "NVj.polarisation",
                "Eion.polarisation",
            ):
                assert terms[key] == 0, key
    return {
        "variant": variant,
        "elliptic_and_X_boundaries": "verified symbolically",
        "reynolds_skew_identity": "verified symbolically",
        "positive_primitives": True,
        "sample_term_magnitudes": magnitudes,
    }
