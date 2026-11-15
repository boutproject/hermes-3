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
    fill_wall_guards,
    reconstruct_reynolds_faces,
)
from energy_analysis.inventory import energy_densities
from energy_analysis.operators import perpendicular_operator
from energy_analysis.transport import (
    face_divergence,
    face_work,
    fv_xz,
    reynolds_face_fluxes,
    weighted_faces,
    xppm_faces,
)


def expressions(variant):
    fields = exact.model(variant)[0]
    phi, psi = fields["phi"], fields["psi"]
    Phi, c, H = phi + psi, s.Integer(2), exact.metric(variant)
    omega = fields["Vort"]
    gradient = H * exact.grad(Phi)
    Gamma = c * s.Matrix(
        [exact.grad(gradient[i]).dot(exact.drift(phi)) for i in range(2)]
    )
    parts = {"omega_vE": omega * exact.drift(phi)}
    parts["R"] = Gamma - parts["omega_vE"]
    rho = 2 * fields["Ni"] + 3 * fields["Nj"]
    ki = fields["NVi"] ** 2 / (4 * fields["Ni"])
    kj = fields["NVj"] ** 2 / (6 * fields["Nj"])
    weights = {
        "Ni": 2 * fields["Ni"] / rho,
        "Nj": s.Rational(3, 2) * fields["Nj"] / rho,
        "NVi": 2 * fields["NVi"] / rho,
        "NVj": s.Rational(3, 2) * fields["NVj"] / rho,
        "E_adv": Phi
        + 2 / rho * (s.Rational(3, 2) * fields["pi"] + ki)
        + s.Rational(3, 2) / rho * (s.Rational(3, 2) * fields["pj"] + kj),
    }
    kinetic = c * exact.grad(Phi).dot(H * exact.grad(Phi)) / 2
    continuum = {k: exact.div(v) for k, v in parts.items()}
    continuum.update(
        DivJReynolds=-exact.div(Gamma),
        DivEExB_perpendicular=exact.div(kinetic * exact.drift(phi)),
        work=Gamma.dot(exact.grad(Phi)),
    )
    for name, weight in weights.items():
        continuum["transport_" + name] = -exact.div(weight * Gamma)
    return fields, weights, continuum


def sample(expression, x, z, y, time):
    value = s.lambdify(exact.COORDINATES, expression, "numpy", cse=True)(
        x[:, None], y, z[None, :], time
    )
    return xr.DataArray(
        np.broadcast_to(value, (len(x), len(z))).copy(), dims=("x", "z")
    )


def discrete(fields, weights, wall_weights, geometry):
    wall = SimpleNamespace(mxg=2)
    phi = fill_wall_guards(fields["phi"], wall, "zero")
    psi = fill_wall_guards(fields["psi"], wall, "pressure")
    omega = fill_wall_guards(fields["Vort"], wall, "neumann")
    c, Phi = 2.0, phi + psi
    faces = reynolds_face_fluxes(c, omega, phi, psi, geometry)
    advective = xppm_faces(omega, phi, geometry)
    parts = {"omega_vE": advective, "R": {d: faces[d] - advective[d] for d in faces}}
    kinetic = perpendicular_operator(c, Phi, geometry).kperp
    result = {
        **{k: face_divergence(v, geometry) for k, v in parts.items()},
    }
    result.update(
        DivJReynolds=-face_divergence(faces, geometry),
        DivEExB_perpendicular=face_divergence(
            xppm_faces(kinetic, phi, geometry), geometry
        ),
        work=face_work(faces, Phi, geometry),
    )
    for name, weight in weights.items():
        transported = weighted_faces(faces, weight, wall_weight=wall_weights[name])
        result["transport_" + name] = -face_divergence(transported, geometry)
    return result


def errors(actual, expected, n, x, z):
    error = np.asarray(actual - expected)[2:-2]
    maximum = np.unravel_index(np.argmax(abs(error)), error.shape)
    return {
        "all_l2": float(np.sqrt(np.mean(error**2))),
        "wall_l2": float(np.sqrt(np.mean(error[[0, -1]] ** 2))),
        "wall_maximum": float(abs(error[[0, -1]]).max()),
        "outer_three_l2": float(
            np.sqrt(np.mean(np.concatenate((error[:3], error[-3:])) ** 2))
        ),
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
                + fill_wall_guards(fields["psi"], wall, "pressure")
                + 2
                / rho
                * (
                    1.5 * guarded["Ni"] * temperature
                    + guarded["NVi"] ** 2 / (4 * guarded["Ni"])
                )
                + 1.5
                / rho
                * (
                    1.5 * guarded["Nj"] * temperature
                    + guarded["NVj"] ** 2 / (6 * guarded["Nj"])
                )
            )
            wall_weights = {name: xr.zeros_like(temperature) for name in w}
            Phi = fill_wall_guards(fields["phi"], wall, "zero") + fill_wall_guards(
                fields["psi"], wall, "pressure"
            )
            for face, cell in ((1, 2), (n + 1, n + 1)):

                def midpoint(value, face=face):
                    return (
                        value.isel(x=face, drop=True)
                        + value.isel(x=face + 1, drop=True)
                    ) / 2

                thermal = midpoint(temperature).mean("z")
                energy = midpoint(Phi).mean("z")
                for species, mass, charge in (("i", 2, 1), ("j", 3, 2)):
                    density = guarded["N" + species]
                    partition = (mass * density / rho).isel(x=cell, drop=True).mean(
                        "z"
                    ) / charge
                    velocity = midpoint(
                        guarded["NV" + species] / (mass * density)
                    ).mean("z")
                    wall_weights["N" + species][{"x": face}] = partition
                    wall_weights["NV" + species][{"x": face}] = (
                        mass * partition * velocity
                    )
                    energy += partition * (1.5 * thermal + 0.5 * mass * velocity**2)
                wall_weights["E_adv"][{"x": face}] = energy
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
            for reconstruction in ("mixed",):
                actual = discrete(fields, w, wall_weights, g)
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
                for region in ("all", "interior", "wall", "outer_three"):
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
        faces = reconstruct_reynolds_faces(run, ds.phi, psi)
        work = face_work(faces, ds.phi + psi, ds)
        energy = face_divergence(
            xppm_faces(ds.Eion, ds.phi, ds, mxg=run.mxg),
            ds,
        )
        chi_faces = fv_xz(
            c,
            ds.ion_chi,
            ds,
            mxg=run.mxg,
            boundary_flux=True,
        )
        polarisation = {
            "polarisation_flux_" + d: -faces[d] - chi_faces["flux_" + d]
            for d in ("x", "z")
        }
        agreement = compare_flux_diagnostics(
            run,
            xr.Dataset(
                {
                    "DivJReynolds": -face_divergence(faces, ds),
                    "DivEExB": energy,
                    "reynolds_work": work,
                    **polarisation,
                }
            ),
        )
        required = (
            "DivJReynolds",
            "DivEExB",
            "reynolds_work",
            "polarisation_flux_x",
            "polarisation_flux_z",
        )
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
