"""Fixed physical poloidal-band wall and interior monitors from native dumps."""

import json
from contextlib import ExitStack
from pathlib import Path

import numpy as np

from .dataset import AuditError

BANDS = {"16-31": (16, 31), "32-63": (32, 63)}


def band_rms(field, low, high):
    """Parseval RMS in inclusive mode numbers, with Z as the final axis."""
    nz = field.shape[-1]
    if not 0 <= low <= high <= nz // 2:
        raise AuditError("Requested physical band is not resolved on this grid")
    transform = np.fft.rfft(field, axis=-1) / nz
    weights = np.full(nz // 2 + 1, 2.0)
    weights[0] = 1.0
    if nz % 2 == 0:
        weights[-1] = 1.0
    return np.sqrt((weights * abs(transform) ** 2)[..., low : high + 1].sum(axis=-1))


def e_fold_width(profile, dx_mm, interior):
    """First inward log-interpolated RMS/e crossing, with wall dominance flag."""
    result = {"mm": None, "cells": None, "wall_dominated": bool(profile[0] > interior)}
    if profile[0] <= 0:
        return result
    crossings = np.flatnonzero(profile[1 : len(profile) // 2] <= profile[0] / np.e) + 1
    if len(crossings):
        j = int(crossings[0])
        left, right = np.log(np.maximum(profile[[j - 1, j]], np.finfo(float).tiny))
        cells = float(j - 0.5 + (np.log(profile[0]) - 1 - left) / (right - left))
        result.update(cells=cells, mm=cells * dx_mm)
    return result


def analyse_native_bands(
    folder,
    *,
    radial_length_mm=50.0,
    poloidal_length_m=0.05,
    interior_exclusion_mm=3.125,
):
    """Stream one X-decomposed drift-plane output at a time, retaining profiles.

    Each rank owns complete Z lines and singleton Y. Lengths describe the
    physical domain and must match the input. Amplitudes use first-cell centres;
    e-fold distances include their half-cell offset from the physical wall.
    Wall/interior ratios are deliberately absent from acceptance criteria.
    """
    from netCDF4 import Dataset

    folder = Path(folder)
    paths = sorted(
        folder.glob("BOUT.dmp.*.nc"), key=lambda p: int(p.name.split(".")[2])
    )
    if not paths:
        raise AuditError(f"No native BOUT++ dumps in {folder}")
    report = {
        "bands": BANDS,
        "radial_length_mm": radial_length_mm,
        "poloidal_length_m": poloidal_length_m,
        "interior_exclusion_mm": interior_exclusion_mm,
        "amplitude_units": "native normalised units",
        "amplitude_location": "first physical cell centres",
    }
    with ExitStack() as stack:
        dumps = [stack.enter_context(Dataset(path)) for path in paths]
        times = np.asarray(dumps[0]["t_array"][:])
        for dump in dumps:
            if int(dump["evolve_vorticity_ion_energy_version"][...]) not in (
                4,
                5,
                6,
                7,
            ):
                raise AuditError(
                    "Fixed-band monitor requires ion-energy version 4, 5, 6 or 7"
                )
            if dump["Vort"].shape[2] != 1 or int(dump["NYPE"][...]) != 1:
                raise AuditError("Band monitor requires a singleton-Y drift plane")
            if int(dump["NXPE"][...]) != len(dumps):
                raise AuditError("Missing native rank dumps")
            np.testing.assert_array_equal(times, dump["t_array"][:])
        guards = int(dumps[0]["MXG"][...])
        nx = sum(dump["Vort"].shape[1] - 2 * guards for dump in dumps)
        dx = radial_length_mm / nx
        cut = round(interior_exclusion_mm / dx)
        if cut < 1 or 2 * cut >= nx:
            raise AuditError("Interior exclusion leaves no resolved interior")
        report.update(
            endpoint=float(times[-1]),
            n=nx,
            cell_width_mm=dx,
            time=times.tolist(),
            distance_mm=((np.arange(nx) + 0.5) * dx).tolist(),
        )
        for variable, names in {"vorticity": ("Vort",), "Pi": ("Nh+", "Tion")}.items():
            data = {
                band: {
                    key: []
                    for key in (
                        "inner_wall_RMS",
                        "outer_wall_RMS",
                        "combined_wall_RMS",
                        "interior_RMS",
                        "profile",
                        "inner_width",
                        "outer_width",
                    )
                }
                for band in BANDS
            }
            for j in range(len(times)):
                fields = []
                for dump in dumps:
                    field = np.asarray(dump[names[0]][j, guards:-guards, 0, :])
                    if len(names) == 2:
                        field *= np.asarray(dump[names[1]][j, guards:-guards, 0, :])
                    if not np.isfinite(field).all():
                        raise AuditError("Nonfinite or incomplete native output")
                    fields.append(field)
                field = np.concatenate(fields, axis=0)
                for band, (low, high) in BANDS.items():
                    profile = band_rms(field, low, high)
                    interior = float(np.sqrt(np.mean(profile[cut:-cut] ** 2)))
                    row = data[band]
                    for key, value in (
                        ("inner_wall_RMS", float(profile[0])),
                        ("outer_wall_RMS", float(profile[-1])),
                        (
                            "combined_wall_RMS",
                            float(np.sqrt(np.mean(profile[[0, -1]] ** 2))),
                        ),
                        ("interior_RMS", interior),
                        ("profile", profile.tolist()),
                        ("inner_width", e_fold_width(profile, dx, interior)),
                        ("outer_width", e_fold_width(profile[::-1], dx, interior)),
                    ):
                        row[key].append(value)
            for band, row in data.items():
                row["maximum_individual_wall_RMS"] = max(
                    max(row["inner_wall_RMS"]), max(row["outer_wall_RMS"])
                )
                row["wavenumber_cycles_per_m"] = [
                    mode / poloidal_length_m for mode in BANDS[band]
                ]
            report[variable] = data
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--radial-length-mm", type=float, default=50.0)
    parser.add_argument("--poloidal-length-m", type=float, default=0.05)
    parser.add_argument("--interior-exclusion-mm", type=float, default=3.125)
    args = parser.parse_args()
    report = analyse_native_bands(
        args.folder,
        radial_length_mm=args.radial_length_mm,
        poloidal_length_m=args.poloidal_length_m,
        interior_exclusion_mm=args.interior_exclusion_mm,
    )
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
