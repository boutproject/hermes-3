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


def metric(variant="parallel_current_wall"):
    sigma = s.Rational(1, 4)
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
    return 1 - s.cos(2 * a)


@lru_cache(maxsize=len(VARIANTS))
def model(variant="parallel_current_wall"):
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
    nu = r(2, 100)
    # Wall-uniform density and velocity isolate the parallel-current operator
    # comparison; temperature retains the local pressure-potential fluctuation.
    wall_envelope = s.sin(a) ** 2
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
    amplitude = r(2, 10)
    sigma = H[0, 1]
    wave_x, wave_z = s.pi / LENGTH, 2 * s.pi / LENGTH
    sheared_phase = theta - sigma * wave_z * X
    # L(Phi) is proportional to sin(a)^3: its value and normal derivative
    # vanish at both fixed walls, matching the closed viscous-current faces.
    envelope = amplitude * (s.sin(a) - r(5, 39) * s.sin(3 * a))
    fluctuation = envelope * s.cos(b + tau) * s.cos(sheared_phase)
    amp_psi = r(2, 100)
    shape = pressure_shape(a, variant)
    pressure_fluctuation = amp_psi * shape * s.cos(sheared_phase)
    psi = r(12, 10) + r(4, 100) * s.cos(b - tau) + pressure_fluctuation
    phi = fluctuation - pressure_fluctuation
    ti = common_temperature(psi, ni, nj)
    pi, pj = ni * ti, nj * ti
    ne = ni + 2 * nj
    ji = ni * vi + 2 * nj * vj
    curvature = r(1, 10)
    # A finite first-Y-harmonic total pressure permits an analytic periodic
    # current while keeping magnetic curvature and electron work nonzero.
    total_pressure = 4 + r(1, 50) * s.cos(2 * a) * s.sin(theta) * s.cos(b - tau)
    pe = total_pressure - pi - pj
    mi, mj = 2 * ni * vi, 3 * nj * vj
    ki, kj = ni * vi**2, r(3, 2) * nj * vj**2
    rho, pressure, Phi, c = 2 * ni + 3 * nj, pi + pj, phi + psi, s.Integer(2)
    omega = c * lap(Phi, H)
    kperp = c * grad(Phi).dot(H * grad(Phi)) / 2
    energy = r(3, 2) * pressure + ki + kj + kperp
    # The full current is -c*advection(H*grad(Phi), phi). Separate its quadratic
    # Y dependence from the first harmonic to integrate Jparallel analytically.
    reynolds_self = -c * s.Matrix(
        [advection(q, fluctuation) for q in H * grad(fluctuation)]
    )
    reynolds_cross = c * s.Matrix(
        [advection(q, pressure_fluctuation) for q in H * grad(fluctuation)]
    )
    reynolds = reynolds_self + reynolds_cross
    chi_harmonic = r(1, 100) * s.cos(2 * a) * s.cos(b - tau) * s.cos(theta) + r(
        1, 200
    ) * s.cos(a) * s.sin(b + tau) * s.cos(2 * theta)
    # Solve the Y-mean continuum Reynolds equation in the sheared coordinates.
    # Q is its finite Fourier particular solution. The homogeneous term imposes
    # coordinate-Neumann AC guards; all terms have zero Z mean at the walls.
    Q = amplitude**2 * wave_x**2 * (r(2, 39) * s.sin(2 * a) - r(5, 312) * s.sin(4 * a))
    wall_derivative = amplitude**2 * wave_x**3 / 26
    harmonic = 2 * wave_z
    centred = harmonic * (X - LENGTH / 2)
    half_width = harmonic * LENGTH / 2
    normalisation = s.exp(half_width) + s.exp(-half_width)
    even = (s.exp(centred) + s.exp(-centred)) / normalisation
    odd = (s.exp(centred) - s.exp(-centred)) / normalisation
    denominator = harmonic * (1 + sigma**2)
    chi_mean = Q * s.sin(2 * sheared_phase) - wall_derivative / denominator * (
        sigma * even * s.cos(2 * sheared_phase) + odd * s.sin(2 * sheared_phase)
    )
    chi = chi_harmonic + chi_mean
    mean_reynolds = (
        c
        * wave_z
        / 4
        * (
            envelope * s.diff(envelope, X, 3)
            - s.diff(envelope, X) * s.diff(envelope, X, 2)
        )
        * s.sin(2 * sheared_phase)
    )
    current_derivative = (
        c * lap(chi_harmonic, H)
        - div(reynolds_cross)
        - curvature * s.diff(total_pressure, Z)
    )
    # The first harmonic is C*cos(b)+S*sin(b). The quadratic part contributes
    # -mean_reynolds*cos(2*(b+tau)); both have periodic elementary integrals.
    cosine = current_derivative.subs(Y, 0)
    sine = current_derivative.subs(Y, LENGTH / 4)
    parallel_current = (cosine * s.sin(b) - sine * s.cos(b)) / wave_z
    parallel_current -= mean_reynolds * s.sin(2 * (b + tau)) / (2 * wave_z)
    je = parallel_current - ji
    ve = -je / ne
    jnu = nu * H * grad(omega)
    jpol = -c * H * grad(chi) + reynolds
    terms = {"Vort.reynolds": div(reynolds)}
    terms["Vort.parallel_current"] = s.diff(parallel_current, Y)
    # Continuum current identity; viscous stress is excluded from chi.
    rhs = {"Vort": c * lap(chi, H)}
    terms["Vort.curvature"] = curvature * s.diff(total_pressure, Z)
    terms["Vort.viscosity"] = div(jnu)
    rhs["Vort"] += terms["Vort.viscosity"]

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
            "Eion.reynolds": -advection(kperp, phi),
            "Eion.parallel_electric": -s.diff(phi * ji, Y),
            "Eion.parallel": -s.diff(
                (r(5, 2) * pi + ki) * vi + (r(5, 2) * pj + kj) * vj, Y
            ),
            # Complete polarisation flux, including finite current divergence.
            "Eion.polarisation": polarisation(
                Phi
                + 2 / rho * (r(3, 2) * pi + ki)
                + r(3, 2) / rho * (r(3, 2) * pj + kj)
            ),
            "Eion.electron_work": -phi * (s.diff(je, Y) + curvature * s.diff(pe, Z)),
        }
    )
    terms["Pe.curvature_work"] = -r(2, 3) * curvature * pe * s.diff(phi, Z)
    terms["Eion.curvature_electric"] = -curvature * s.diff(phi * pressure, Z)
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
        Ve=ve,
        Jpar=parallel_current,
        DivJpar=s.diff(parallel_current, Y),
    )
    sources = {name: s.diff(fields[name], t) - rhs[name] for name in evolved(variant)}
    return fields, rhs, sources, terms


@cache
def evaluator(name, variant="parallel_current_wall"):
    return s.lambdify(COORDINATES, model(variant)[0][name], "numpy", cse=True)


def evaluate(name, x, y, z, time, variant="parallel_current_wall"):
    """Evaluate on broadcastable physical coordinates, including guards."""
    result = evaluator(name, variant)(x, y, z, time)
    return np.broadcast_to(
        result,
        np.broadcast_shapes(np.shape(x), np.shape(y), np.shape(z), np.shape(time)),
    )


def zero(expression):
    """Resolve finite trigonometric identities using independent differentiation."""
    return s.simplify(s.expand(s.expand_trig(expression).rewrite(s.exp))) == 0


def preflight(variant="parallel_current_wall"):
    """Verify recovery, physical current continuity, boundaries and term coverage."""
    fields, rhs, sources, terms = model(variant)
    H, c = metric(variant), s.Integer(2)
    Phi = fields["phi"] + fields["psi"]
    assert zero(fields["Vort"] - c * lap(Phi, H))
    assert zero(
        c * lap(fields["chi"], H)
        - sum(
            terms[name]
            for name in ("Vort.reynolds", "Vort.parallel_current", "Vort.curvature")
        )
    ), "Auxiliary current identity"
    assert zero(rhs["Vort"] - c * lap(fields["chi"], H) - terms["Vort.viscosity"]), (
        "Viscosity must be excluded from chi"
    )
    current = fields["NVi"] / 2 + 2 * fields["NVj"] / 3 - 1836 * fields["NVe"]
    assert zero(current - fields["Jpar"])
    assert zero(s.diff(current, Y) - fields["DivJpar"])
    for face in (0, LENGTH):
        assert zero(fields["phi"].subs(X, face))
        assert not zero(s.diff(fields["phi"], X).subs(X, face))
        assert zero(s.diff(fields["Vort"], X).subs(X, face))
        assert zero((H * grad(fields["Vort"]))[0].subs(X, face)), (
            "Fixed walls close the normal viscous current"
        )
        for name in ("Ni", "Nj", "NVi", "NVj", "Pe", "psi", "chi"):
            assert zero(s.diff(fields[name], X).subs(X, face)), name
        assert not zero(s.diff(fields["psi"], X, 2).subs(X, face))
        assert zero(s.integrate(fields["chi"].subs(X, face), (Z, 0, LENGTH)))
    # psi lies in [1.10,1.30]. Pi <= 4*psi/3, hence Pe >= 4-.02-4*1.30/3.
    assert s.Rational(4) - s.Rational(1, 50) - s.Rational(4, 3) * s.Rational(13, 10) > 0
    magnitudes = {
        key: abs(
            float(
                s.lambdify(COORDINATES, value, "numpy", cse=True)(
                    1.37, 2.91, 4.23, 0.73
                )
            )
        )
        for key, value in terms.items()
    }
    assert all(np.isfinite(value) and value > 1e-12 for value in magnitudes.values()), (
        magnitudes
    )
    assert set(evolved(variant)) <= sources.keys() and set(RECOVERED) <= fields.keys()
    return {
        "elliptic_and_X_boundaries": "verified symbolically",
        "positive_primitives": True,
        "g13": 0.25,
        "viscosity": 0.02,
        "magnetic_curvature": 0.1,
        "sample_term_magnitudes": magnitudes,
    }
