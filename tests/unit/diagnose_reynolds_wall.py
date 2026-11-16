#!/usr/bin/env python3
"""Bounded, solver-free wall-current probes; no production boundary changes.

Use the retained flat-wall state and the independent covariant reference.
Writes area-weighted midpoint-face errors, split divergence errors, radial profiles and
Gauss balances. All-cell error is mean absolute; boundary is max, outer-three
and interior are RMS. Interior excludes three physical cells at either wall.
"""

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import sympy as s
import xarray as xr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tests/python"))
from energy_analysis.transport import (
    face_divergence,
    fv_xz,
    mc_edges,
    neighbor,
    reynolds_face_fluxes,
)
from test_reynolds_transport import reference, wall_curvature_state


def norms(error, *, faces=False):
    a = np.asarray(error)
    a = a[1:-2] if faces else a[2:-2]
    return {
        "all": float(np.mean(abs(a))),
        "boundary": float(np.max(abs(a[[0, -1]]))),
        "outer_three": float(np.sqrt(np.mean(a[[0, 1, 2, -3, -2, -1]] ** 2))),
        "interior": float(np.sqrt(np.mean(a[3:-3] ** 2))),
        "max": float(np.max(abs(a))),
        "max_location": [int(i) for i in np.unravel_index(abs(a).argmax(), a.shape)],
    }


def probe(nx, nz, kind):
    analytic = kind in ("sampled", "continuation", "analytic_faces")
    phi, psi, omega, g, expected, *_ = wall_curvature_state(
        nx, nz, recovered=not analytic, analytic_guards=kind == "continuation"
    )
    h, dz = 1 / nx, 2 * np.pi / nz
    x = (np.arange(nx + 4) - 1.5) * h
    z = np.arange(nz) * dz
    v = reference.reference("wall")

    def sample(expr, dim=None):
        values = s.lambdify((reference.x, reference.z), expr, "numpy", cse=True)(
            (x + (h / 2 if dim == "x" else 0))[:, None],
            (z + (dz / 2 if dim == "z" else 0))[None, :],
        )
        return xr.DataArray(np.broadcast_to(values, phi.shape).copy(), dims=phi.dims)

    # Independent continuum split; sum error fields before evaluating a norm.
    vel = {"x": s.diff(v["phi"], reference.z), "z": -s.diff(v["phi"], reference.x)}
    # Midpoint current times face area, not an exact face integral.
    exact_faces = {
        d: sample(v["expected_" + d], d) * (dz if d == "x" else h) for d in vel
    }
    exact_adv = {
        d: sample(v["omega"] * vel[d], d) * (dz if d == "x" else h) for d in vel
    }
    div_adv = sum(
        s.diff(v["omega"] * vel[d], q) for d, q in zip(vel, (reference.x, reference.z))
    )
    adv_expected = sample(div_adv)
    r_expected = expected - adv_expected
    faces = reynolds_face_fluxes(
        1, omega, phi, psi, g, limiter="central" if kind == "central" else "mc"
    )
    # Resolve the numerical Omega*v_E and R contributions at the same faces.
    xp, xm = lambda f: neighbor(f, "x"), lambda f: neighbor(f, "x", -1)
    zp, zm = lambda f: neighbor(f, "z"), lambda f: neighbor(f, "z", -1)
    corner = (phi + xp(phi) + zp(phi) + xp(zp(phi))) / 4
    velocity = {"x": (corner - zm(corner)) / dz, "z": (xm(corner) - corner) / h}
    adv = {}
    for d, shift in (("x", xp), ("z", zp)):
        left, right = mc_edges(omega, d)
        donor = (
            (omega + shift(omega)) / 2
            if kind == "central"
            else xr.where(velocity[d] > 0, right, shift(left))
        )
        adv[d] = velocity[d] * donor * (dz if d == "x" else h)
    if kind == "analytic_faces":
        faces, adv = exact_faces, exact_adv
    remainder = {d: faces[d] - adv[d] for d in faces}
    divergence = face_divergence(faces, g)
    error = divergence - expected
    adv_error = face_divergence(adv, g) - adv_expected
    r_error = face_divergence(remainder, g) - r_expected
    volume = h * dz
    gauss = float(divergence.isel(x=slice(2, -2)).sum()) * volume
    wall_flux = float((faces["x"].isel(x=-3) - faces["x"].isel(x=1)).sum())
    recovery = omega - fv_xz(1, phi + psi, g).divergence
    # Fourier amplitudes of recovered Phi give a known radial interpolant:
    # cubic splines on physical samples, evaluated with analytic Z dependence.
    # This derivative probe excludes three wall cells; continuation is a diagnostic.
    from scipy.interpolate import CubicSpline

    physical = np.asarray(phi + psi)[2:-2]
    fc = 2 * np.mean(physical * np.cos(z), axis=1)
    fs = 2 * np.mean(physical * np.sin(z), axis=1)
    psi_f = 2 * np.mean(np.asarray(psi)[2:-2] * np.sin(z), axis=1)
    csp = CubicSpline(x[2:-2], fc)
    ssp = CubicSpline(x[2:-2], fs)
    psp = CubicSpline(x[2:-2], psi_f)
    # At each physical x: phi = C cos(z) + (S-Fpsi) sin(z).
    zz = z[None, :]
    C, S, F = [fun(x[2:-2])[:, None] for fun in (csp, ssp, psp)]
    C1, S1, F1 = [fun(x[2:-2], 1)[:, None] for fun in (csp, ssp, psp)]
    C2, S2, F2 = [fun(x[2:-2], 2)[:, None] for fun in (csp, ssp, psp)]
    C3, S3 = [fun(x[2:-2], 3)[:, None] for fun in (csp, ssp)]
    phix = C1 * np.cos(zz) + (S1 - F1) * np.sin(zz)
    phiz = -C * np.sin(zz) + (S - F) * np.cos(zz)
    phixz = -C1 * np.sin(zz) + (S1 - F1) * np.cos(zz)
    phizz = -C * np.cos(zz) - (S - F) * np.sin(zz)
    Pxz = -C1 * np.sin(zz) + S1 * np.cos(zz)
    Pxx = C2 * np.cos(zz) + S2 * np.sin(zz)
    Pzz = -C * np.cos(zz) - S * np.sin(zz)
    Pxxx = C3 * np.cos(zz) + S3 * np.sin(zz)
    Pxzz = -C1 * np.cos(zz) - S1 * np.sin(zz)
    Pxxz = -C2 * np.sin(zz) + S2 * np.cos(zz)
    Pzzz = C * np.sin(zz) - S * np.cos(zz)
    state_div = (
        phixz * Pxx
        + phiz * Pxxx
        - (C2 * np.cos(zz) + (S2 - F2) * np.sin(zz)) * Pxz
        - phix * Pxxz
    )
    state_div += phizz * Pxz + phiz * Pxzz - phixz * Pzz - phix * Pzzz
    record = {
        "nx": nx,
        "nz": nz,
        "kind": kind,
        "divergence": norms(error),
        "advection_divergence": norms(adv_error),
        "R_divergence": norms(r_error),
        "recovery_mismatch": norms(recovery),
        "gauss_residual": gauss - wall_flux,
        "wall_current_max": float(abs(faces["x"].isel(x=[1, -3]) / dz).max()),
        "radial_error_mean_abs": np.mean(abs(np.asarray(error)[2:-2]), axis=1).tolist(),
        "radial_error_max": np.max(abs(np.asarray(error)[2:-2]), axis=1).tolist(),
        "state_interpolant_divergence_interior_rms": float(
            np.sqrt(np.mean((state_div - np.asarray(expected)[2:-2])[3:-3] ** 2))
        ),
        "faces": {},
    }
    for d, area in (("x", dz), ("z", h)):
        face_error = (faces[d] - exact_faces[d]) / area
        record["faces"][d] = {
            "total": norms(face_error, faces=d == "x"),
            "advection": norms((adv[d] - exact_adv[d]) / area, faces=d == "x"),
            "R": norms(
                (remainder[d] - (exact_faces[d] - exact_adv[d])) / area, faces=d == "x"
            ),
        }
    first = np.asarray((faces["x"] - exact_faces["x"]) / dz)[1:4]
    record["first_x_face_errors_over_h"] = (first / h).tolist()
    record["first_cell_divergence_error"] = np.asarray(error)[2].tolist()
    return record, (phi, psi, omega, faces, divergence)


def face_quadrature_probe(n):
    """Separate midpoint/point references from exact face/cell integrals."""
    v = reference.reference("wall")
    x, z = reference.x, reference.z
    h, dz = 1 / n, 2 * np.pi / n
    xc = (np.arange(n) + 0.5) * h
    zc = np.arange(n) * dz

    def values(expr, xx, zz):
        return np.broadcast_to(
            s.lambdify((x, z), expr, "numpy", cse=True)(xx[:, None], zz[None, :]),
            (len(xx), len(zz)),
        )

    px = s.integrate(v["expected_x"], z)
    pz = s.integrate(v["expected_z"], x)
    pc = s.integrate(s.integrate(v["expected"], x), z)
    xfaces = values(px, np.arange(n + 1) * h, zc + dz / 2) - values(
        px, np.arange(n + 1) * h, zc - dz / 2
    )
    zfaces = values(pz, xc + h / 2, zc + dz / 2) - values(pz, xc - h / 2, zc + dz / 2)
    divergence = (np.diff(xfaces, axis=0) + zfaces - np.roll(zfaces, 1, axis=1)) / (
        h * dz
    )
    average = (
        values(pc, xc + h / 2, zc + dz / 2)
        - values(pc, xc - h / 2, zc + dz / 2)
        - values(pc, xc + h / 2, zc - dz / 2)
        + values(pc, xc - h / 2, zc - dz / 2)
    ) / (h * dz)
    point = values(v["expected"], xc, zc)
    residual = float(abs(divergence - average).max())
    if residual > 1e-9:
        raise AssertionError((n, "analytic face/cell integral identity", residual))
    return {
        "n": n,
        "face_integral_vs_cell_average_max": residual,
        "cell_average_vs_point_rms": float(np.sqrt(np.mean((average - point) ** 2))),
        "cell_average_vs_point_boundary_max": float(
            abs(average - point)[[0, -1]].max()
        ),
        "arithmetic_limit": 1e-9,
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    cli.add_argument("--face-quadrature-only", action="store_true")
    args = cli.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    if args.face_quadrature_only:
        probes = [face_quadrature_probe(n) for n in (32, 64, 128)]
        (args.output / "face-quadrature.json").write_text(
            json.dumps(
                {
                    "purpose": "Exact analytic face integrals versus cell averages and point samples",
                    "command": sys.argv,
                    "decomposition": "solver-free, one Python process",
                    "source_sha256": hashlib.sha256(
                        Path(__file__).read_bytes()
                    ).hexdigest(),
                    "elapsed_seconds": time.monotonic() - start,
                    "probes": probes,
                },
                indent=2,
                allow_nan=False,
            )
            + "\n"
        )
        print(probes)
        return
    results = []
    for kind in ("recovered", "sampled", "continuation", "analytic_faces", "central"):
        for direction in ("xz", "x", "z"):
            for n in (32, 64, 128):
                nx, nz = (
                    (n, n)
                    if direction == "xz"
                    else (n, 128)
                    if direction == "x"
                    else (128, n)
                )
                record, fields = probe(nx, nz, kind)
                record["refinement"] = direction
                results.append(record)
                if direction == "xz" and kind in ("recovered", "continuation"):
                    phi, psi, omega, faces, divergence = fields
                    rows = []
                    for ix in range(nx + 4):
                        for iz in range(nz):
                            rows.append(
                                [
                                    ix,
                                    iz,
                                    float(phi[ix, iz]),
                                    float(psi[ix, iz]),
                                    float(omega[ix, iz]),
                                    float(faces["x"][ix, iz]),
                                    float(faces["z"][ix, iz]),
                                    float(divergence[ix, iz]),
                                ]
                            )
                    np.savetxt(
                        args.output / f"cpp-reference-{kind}-{n}.csv",
                        rows,
                        delimiter=",",
                        fmt="%.17g",
                    )
        print(
            kind,
            [
                (r["nx"], r["divergence"]["boundary"], r["divergence"]["all"])
                for r in results
                if r["kind"] == kind and r["refinement"] == "xz"
            ],
            flush=True,
        )
    manifest = {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (
            Path(__file__),
            ROOT / "tests/python/test_reynolds_transport.py",
            ROOT / "tests/unit/generate_reynolds_reference.py",
            ROOT / "tools/energy_analysis/transport.py",
        )
    }
    (args.output / "matrix.json").write_text(
        json.dumps(
            {
                "purpose": __doc__,
                "command": sys.argv,
                "decomposition": "solver-free, one Python process",
                "elapsed_seconds": time.monotonic() - start,
                "source_sha256": manifest,
                "probes": results,
            },
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
