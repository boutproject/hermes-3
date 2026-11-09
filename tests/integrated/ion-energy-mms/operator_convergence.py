#!/usr/bin/env python3
"""Attribute spatial MMS errors on analytic X-Z planes; no Hermes run required.

Sample Y and time explicitly, retaining analytic Z periodicity and the fixed
MMS wall guards. --run additionally checks reconstructed C++ diagnostics.
This diagnostic does not alter the manufactured sources or convergence gates.
"""

import argparse
import json
import sys
from itertools import pairwise
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import sympy as s
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
import manufactured as exact
from energy_analysis import open_run
from energy_analysis.fluxes import (
    compare_flux_diagnostics,
    fill_transport_guards,
    fill_wall_guards,
)
from energy_analysis.inventory import energy_densities
from energy_analysis.operators import perpendicular_operator
from energy_analysis.transport import (
    arakawa,
    central_advection_faces,
    face_divergence,
    face_work,
    fv_xz,
    neighbor,
    reynolds_face_fluxes,
    reynolds_lambda,
    solenoidal_reynolds_faces,
    xppm_faces,
)


def expressions(variant):
    fields = exact.model(variant)[0]
    phi, psi = fields["phi"], fields["psi"]
    Phi, c, H = phi + psi, s.Integer(2), exact.metric(variant)
    D = H.det()
    lam = c * D * exact.div(H * exact.grad(phi) / D) / 2
    Q = exact.advection(psi, phi)
    omega = fields["Vort"]
    split = {
        "omega_vE": omega * exact.drift(phi) / 2,
        "Q": c * H * exact.grad(Q) / 2,
        "Lambda_U": lam * exact.drift(Phi),
    }
    stream = -c * exact.grad(phi).dot(H * exact.grad(Phi)) / 2
    curl = exact.drift(stream)
    Gamma = sum(split.values(), s.zeros(2, 1)) + curl
    rho = 2 * fields["Ni"] + 3 * fields["Nj"]
    ki = fields["NVi"] ** 2 / (4 * fields["Ni"])
    kj = fields["NVj"] ** 2 / (6 * fields["Nj"])
    weights = {
        "Ni": 2 * fields["Ni"] / rho,
        "Nj": s.Rational(3, 2) * fields["Nj"] / rho,
        "NVi": 2 * fields["NVi"] / rho,
        "NVj": s.Rational(3, 2) * fields["NVj"] / rho,
        "E_adv": phi
        + 2 / rho * (s.Rational(5, 2) * fields["pi"] + ki)
        + s.Rational(3, 2) / rho * (s.Rational(5, 2) * fields["pj"] + kj),
    }
    kinetic = c * exact.grad(Phi).dot(H * exact.grad(Phi)) / 2
    continuum = {"Lambda": lam, **{k: exact.div(v) for k, v in split.items()}}
    continuum.update(
        DivJReynolds=-exact.div(Gamma),
        DivEReynolds=exact.div(kinetic * exact.drift(phi)),
        work=Gamma.dot(exact.grad(Phi)),
    )
    for name, weight in weights.items():
        continuum["transport_" + name] = -exact.div(weight * Gamma)
        continuum["S_grad_" + name] = curl.dot(exact.grad(weight))
    return fields, weights, continuum


def sample(expression, x, z, y, time):
    value = s.lambdify(exact.COORDINATES, expression, "numpy", cse=True)(
        x[:, None], y, z[None, :], time
    )
    return xr.DataArray(
        np.broadcast_to(value, (len(x), len(z))).copy(), dims=("x", "z")
    )


def discrete(fields, weights, geometry, reconstruction):
    wall = SimpleNamespace(mxg=2)
    phi = fill_wall_guards(fields["phi"], wall, "zero")
    psi = fill_wall_guards(fields["psi"], wall, "pressure")
    omega = fill_wall_guards(fields["Vort"], wall, "neumann")
    c, Phi = 2.0, phi + psi
    lam = reynolds_lambda(c, phi, geometry)
    guarded_lam = fill_wall_guards(lam, wall, "free_o2")
    Q = fill_wall_guards(arakawa(phi, psi, geometry.dx, geometry.dz), wall, "zero")
    limiter = None if reconstruction == "central" else "mc"
    parts = {
        "omega_vE": xppm_faces(omega / 2, phi, geometry, limiter=limiter),
        "Q": {
            d: fv_xz(c / 2, Q, geometry, boundary_flux=False)["flux_" + d]
            for d in ("x", "z")
        },
        "Lambda_U": xppm_faces(guarded_lam, Phi, geometry, limiter=limiter),
    }
    faces = reynolds_face_fluxes(
        c, omega, phi, psi, geometry, reconstruction=reconstruction
    )
    curl = solenoidal_reynolds_faces(c, phi, Phi, geometry, walls=True)
    kinetic = perpendicular_operator(c, Phi, geometry).kperp
    result = {
        "Lambda": lam,
        **{k: face_divergence(v, geometry) for k, v in parts.items()},
    }
    result.update(
        DivJReynolds=-face_divergence(faces, geometry),
        DivEReynolds=face_divergence(
            central_advection_faces(kinetic, phi, geometry), geometry
        ),
        work=face_work(faces, Phi, geometry),
    )
    for name, weight in weights.items():
        weighted = lambda flux, weight=weight: {
            d: f * (weight + neighbor(weight, d)) / 2 for d, f in flux.items()
        }
        result["transport_" + name] = -face_divergence(weighted(faces), geometry)
        result["S_grad_" + name] = face_divergence(weighted(curl), geometry)
    return result


def errors(actual, expected, n, x, z):
    error = np.asarray(actual - expected)[2:-2]
    maximum = np.unravel_index(np.argmax(abs(error)), error.shape)
    return {
        "all_l2": float(np.sqrt(np.mean(error**2))),
        "interior_l2": float(np.sqrt(np.mean(error[2:-2] ** 2))),
        "maximum": float(abs(error[maximum])),
        "maximum_cell": [int(i) for i in maximum],
        "maximum_position": [float(x[maximum[0] + 2]), float(z[maximum[1]])],
        "maximum_at_wall": bool(maximum[0] < 2 or maximum[0] >= n - 2),
    }


def convergence(args):
    report = {"y": args.y, "time": args.time, "sizes": args.sizes, "variants": {}}
    for variant in args.variants:
        analytic, weights, continuum = expressions(variant)
        reference = {}
        for n in args.sizes:
            h = float(exact.LENGTH) / n
            x, z = (np.arange(n + 4) - 1.5) * h, (np.arange(n) + 0.5) * h
            evaluate = lambda f, x=x, z=z: sample(f, x, z, args.y, args.time)
            fields = {name: evaluate(expr) for name, expr in analytic.items()}
            w = {name: evaluate(expr) for name, expr in weights.items()}
            # Density and momentum primitives have radial Neumann conditions.
            # Recover pressure-derived energy weights with the component's
            # separately imposed total-pressure and potential wall guards.
            wall = SimpleNamespace(mxg=2)
            for name, value in w.items():
                if name != "E_adv":
                    w[name] = fill_wall_guards(value, wall, "neumann")
            guarded = {
                name: fill_wall_guards(fields[name], wall, "neumann")
                for name in ("Ni", "Nj", "NVi", "NVj")
            }
            pressure = fill_wall_guards(fields["Pi"], wall, "pressure")
            temperature = pressure / (guarded["Ni"] + guarded["Nj"])
            rho = 2 * guarded["Ni"] + 3 * guarded["Nj"]
            w["E_adv"] = (
                fill_wall_guards(fields["phi"], wall, "zero")
                + 2
                / rho
                * (
                    2.5 * guarded["Ni"] * temperature
                    + guarded["NVi"] ** 2 / (4 * guarded["Ni"])
                )
                + 1.5
                / rho
                * (
                    2.5 * guarded["Nj"] * temperature
                    + guarded["NVj"] ** 2 / (6 * guarded["Nj"])
                )
            )
            H = exact.metric(variant)
            g = xr.Dataset(
                {
                    "J": 1.0,
                    "Bxy": 1.0,
                    "dx": h,
                    "dy": 1.0,
                    "dz": h,
                    "g11": float(H[0, 0]),
                    "g13": float(H[0, 1]),
                    "g33": float(H[1, 1]),
                }
            )
            expected = {name: evaluate(expr) for name, expr in continuum.items()}
            for reconstruction in args.reconstruction:
                actual = discrete(fields, w, g, reconstruction)
                terms = reference.setdefault(reconstruction, {})
                for name in continuum:
                    terms.setdefault(name, {})[str(n)] = errors(
                        actual[name], expected[name], n, x, z
                    )
        report["variants"][variant] = reference
        for reconstruction, terms in reference.items():
            print(
                f"\n{variant}, {reconstruction}: all / interior L2 orders; finest maximum"
            )
            for name, values in terms.items():
                orders = {}
                for region in ("all", "interior"):
                    err = np.array([values[str(n)][region + "_l2"] for n in args.sizes])
                    # Zero-phi variants have identically zero Reynolds terms.
                    # Their orders are undefined, rather than NaN/Infinity.
                    orders[region] = [
                        float(np.log2(a / b)) if a > 0 and b > 0 else None
                        for a, b in pairwise(err)
                    ]
                terms[name] = {"errors": values, "orders": orders}
                finest = values[str(args.sizes[-1])]
                print(
                    f"{name:20s} {np.array(orders['all'])} / {np.array(orders['interior'])}; "
                    f"{finest['maximum']:.4g} at {finest['maximum_position']} "
                    f"({'wall band' if finest['maximum_at_wall'] else 'interior'})"
                )
    if args.run:
        run = open_run(args.run)
        # Compare the X-Z operator only: MMS parallel compression and sound
        # speed are outside the full-RHS audit contract and this attribution.
        ds = run.data
        density = energy_densities(run)
        psi = density.pressure_potential
        c = run.scalar_option("average_atomic_mass") / ds.Bxy**2
        omega = fill_wall_guards(ds.Vort, run, "neumann")
        reconstruction = (
            str(ds.options[run.model]["reynolds_reconstruction"])
            if "reynolds_reconstruction" in ds.options[run.model]
            else "xppm"
        )
        faces = reynolds_face_fluxes(
            c,
            omega,
            ds.phi,
            psi,
            ds,
            mxg=run.mxg,
            reconstruction=reconstruction,
            fill_guards=lambda f: fill_transport_guards(f, run, False, True),
        )
        work = face_work(faces, ds.phi + psi, ds)
        energy = face_divergence(
            central_advection_faces(density.ion_perpendicular, ds.phi, ds, mxg=run.mxg),
            ds,
        )
        agreement = compare_flux_diagnostics(
            run,
            xr.Dataset(
                {
                    "DivJReynolds": -face_divergence(faces, ds),
                    "DivEReynolds": energy,
                    "reynolds_work": work,
                    "reynolds_heating": work - energy,
                }
            ),
        )
        required = ("DivJReynolds", "DivEReynolds", "reynolds_work", "reynolds_heating")
        for name in required:
            if (
                not agreement[name]["available"]
                or agreement[name]["max_absolute_error_normalised"] > 1e-10
            ):
                raise AssertionError(
                    f"C++ mirror disagreement: {name}: {agreement[name]}"
                )
        report["cpp_agreement"] = agreement
    return report


def main():
    cli = argparse.ArgumentParser(__doc__)
    cli.add_argument(
        "--variants",
        nargs="+",
        choices=exact.VARIANTS,
        default=("baseline", "curvature"),
    )
    cli.add_argument("--sizes", nargs="+", type=int, default=(16, 32, 64))
    cli.add_argument(
        "--reconstruction",
        nargs="+",
        choices=("xppm", "central"),
        default=("xppm", "central"),
    )
    cli.add_argument("--y", type=float, default=2.91)
    cli.add_argument("--time", type=float, default=0.1)
    cli.add_argument(
        "--run", type=Path, help="Optional diagnosed MMS output for C++ agreement"
    )
    cli.add_argument(
        "--output", type=Path, help="Save machine-readable errors and orders"
    )
    args = cli.parse_args()
    if (
        len(args.sizes) < 2
        or any(n < 8 for n in args.sizes)
        or any(b != 2 * a for a, b in zip(args.sizes, args.sizes[1:]))
    ):
        cli.error("--sizes requires at least two successive doubled grids, each >= 8")
    report = convergence(args)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
