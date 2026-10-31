"""Independent MMS errors for evolved and recovered fields, including X faces."""

import csv
import json
import sys
from itertools import pairwise
from pathlib import Path

import numpy as np
import xhermes
from manufactured import EVOLVED, LENGTH, RECOVERED, evaluate, metric

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
from energy_analysis.operators import perpendicular_operator
from energy_analysis.options import physics_options_file


def require(condition, message):
    if not condition:
        raise ValueError(message)


def open_case(folder):
    folder = Path(folder)
    options_path = folder / "BOUT.inp"
    if (folder / "BOUT.settings").exists():
        options_path = folder / "BOUT.settings"
    with physics_options_file(options_path) as physics_path:
        ds = xhermes.open_hermesdataset(
            datapath=str(folder / "BOUT.dmp.*.nc"),
            inputfilepath=str(physics_path),
            unnormalise=False,
            keep_xboundaries=True,
            keep_yboundaries=False,
            chunks={"t": 1},
            info=False,
        )
    ds.options.filename = str(options_path)
    opts = ds.options
    require(not opts.get_bool("periodicX"), "MMS requires physical X boundaries")
    require(
        opts["mesh:paralleltransform:type"] == "identity",
        "MMS requires an identity parallel transform",
    )
    require(
        not opts["evolve_vorticity_ion_energy"].get_bool("phi_boundary_relax"),
        "MMS requires fixed phi boundaries",
    )
    variant = opts["mms:variant"] if "mms:variant" in opts else "baseline"  # noqa: SIM401 -- BoutOptions has no get()
    H = np.asarray(metric(variant)).astype(float)
    inverse = np.linalg.inv(H)
    expected_metric = {
        "g11": H[0, 0],
        "g33": H[1, 1],
        "g13": H[0, 1],
        "g_11": inverse[0, 0],
        "g_33": inverse[1, 1],
        "g_13": inverse[0, 1],
        "g22": 1,
        "g_22": 1,
        "g12": 0,
        "g23": 0,
        "g_12": 0,
        "g_23": 0,
        "J": 1,
        "Bxy": 1,
        "bxcvx": 0,
        "bxcvy": 0,
        "bxcvz": 0,
    }
    for key, value in expected_metric.items():
        require(
            np.isclose(
                opts["mesh"].evaluate_scalar(key), value, rtol=1e-12, atol=1e-14
            ),
            f"Metric/curvature {key} differs from analytic {variant} slab",
        )
    for setting in ("symmetricGlobalX", "symmetricGlobalY"):
        require(opts["mesh"].get_bool(setting), f"Require mesh:{setting}=true")
    mxg = int(ds.metadata["MXG"])
    nx = int(opts["mesh"].evaluate_scalar("nx")) - 2 * mxg
    ny, nz = (int(opts["mesh"].evaluate_scalar(k)) for k in ("ny", "nz"))
    require(
        ds.sizes["x"] == nx + 2 * mxg and ds.sizes["y"] == ny and ds.sizes["z"] == nz,
        "Unexpected retained guards, duplicate cells or nonperiodic Y topology",
    )
    require(
        int(ds.metadata["ixseps1"]) >= nx + 2 * mxg
        and int(ds.metadata["ixseps2"]) >= nx + 2 * mxg
        and int(ds.metadata["jyseps1_1"]) == -1
        and int(ds.metadata["jyseps2_2"]) == ny - 1,
        "Y is not a closed periodic slab",
    )
    length = float(LENGTH)
    for name, value in (("dx", length / nx), ("dy", length / ny), ("dz", length / nz)):
        require(
            np.allclose(ds[name].values, value, rtol=1e-12, atol=1e-14),
            f"Inconsistent physical coordinate spacing {name}",
        )
    times = ds.t.values
    expected_step = float(opts["solver"].evaluate_scalar("output_step"))
    nout = int(opts["solver"].evaluate_scalar("nout"))
    require(
        len(times) == nout + 1
        and np.allclose(
            times, np.arange(nout + 1) * expected_step, rtol=1e-10, atol=1e-12
        ),
        "Incomplete or irregular output time axis",
    )
    # FieldFactory uses GlobalX/GlobalY cell centres with symmetricGlobal*=true;
    # Z has integer-index sampling. Construct coordinates from resolved options,
    # independently of xhermes' display coordinate labels.
    x = (np.arange(nx + 2 * mxg) - mxg + 0.5) * length / nx
    y = (np.arange(ny) + 0.5) * length / ny
    z = np.arange(nz) * length / nz
    return ds, mxg, (nx, ny, nz), (x[:, None, None], y[None, :, None], z[None, None, :])


def fields_at(ds, index):
    fields = {name: ds[name].isel(t=index) for name in EVOLVED}
    ti = ds.Tion.isel(t=index)
    ni, nj = fields["Ni"], fields["Nj"]
    fields.update(
        Pi=(ni + nj) * ti,
        pi=ni * ti,
        pj=nj * ti,
        Tion=ti,
        phi=ds.phi.isel(t=index),
        psi=(2 * ni + 1.5 * nj) * ti / 2,
    )
    if "ddt(Phi)" in ds:
        fields["chi"] = ds["ddt(Phi)"].isel(t=index)
    return fields


def norms(error, volume):
    require(np.isfinite(error).all(), "Nonfinite MMS error")
    return float(np.sqrt(np.sum(volume * error**2) / np.sum(volume))), float(
        np.max(np.abs(error))
    )


def analyse_case(folder, label=None):
    ds, mxg, grid, coordinates = open_case(folder)
    variant = ds.options["mms:variant"] if "mms:variant" in ds.options else "baseline"  # noqa: SIM401 -- BoutOptions has no get()

    def reference_field(name, x, y, z, time):
        return evaluate(name, x, y, z, time, variant)

    selection = {"x": slice(mxg, -mxg)}
    volume = (
        (ds.J * ds.dx * ds.dy * ds.dz)
        .broadcast_like(ds.Ni.isel(t=0))
        .transpose("x", "y", "z")
    )
    volume = volume.isel(selection).values
    require(np.isfinite(volume).all() and (volume > 0).all(), "Invalid cell volumes")
    label = label or Path(folder).name
    rows, residuals = [], []
    for index, time in enumerate(ds.t.values):
        fields = fields_at(ds, index)
        for name, field in fields.items():
            numerical = field.transpose("x", "y", "z").values
            exact = reference_field(name, *coordinates, float(time))
            full_error = numerical - exact
            error, reference = full_error[mxg:-mxg], exact[mxg:-mxg]
            l2, linf = norms(error, volume)
            ref_l2, ref_linf = norms(reference, volume)
            native_discrepancy = None
            if name in EVOLVED:
                native = (
                    ds["E_" + name]
                    .isel(t=index, **selection)
                    .transpose("x", "y", "z")
                    .values
                )
                native_discrepancy = float(np.max(np.abs(native - error)))
                require(
                    native_discrepancy < 2e-11,
                    f"Native/independent reference mismatch {name} at t={time}: {native_discrepancy}",
                )
            boundary_value_error, boundary_derivative_error = None, None
            if name in ("Vort", "Eion", "phi", "chi"):
                differences = []
                for guard, cell, face in (
                    (mxg - 1, mxg, 0),
                    (-mxg, -mxg - 1, float(LENGTH)),
                ):
                    target = reference_field(
                        name, face, coordinates[1], coordinates[2], float(time)
                    )
                    differences.append(
                        np.abs((numerical[guard] + numerical[cell]) / 2 - target)
                    )
                boundary_value_error = float(np.max(differences))
            else:
                h = float(LENGTH) / grid[0]
                boundary_derivative_error = float(
                    max(
                        np.max(np.abs(numerical[mxg] - numerical[mxg - 1])) / h,
                        np.max(np.abs(numerical[-mxg] - numerical[-mxg - 1])) / h,
                    )
                )
            if name in ("Vort", "Eion", "phi"):
                require(
                    boundary_value_error < 2e-10,
                    f"Incorrect {name} X face value at t={time}: {boundary_value_error}",
                )
            elif name != "chi":
                require(
                    boundary_derivative_error < 2e-10,
                    f"Incorrect {name} Neumann X guards at t={time}: {boundary_derivative_error}",
                )
            position = np.unravel_index(np.argmax(np.abs(error)), error.shape)
            rows.append(
                {
                    "case": label,
                    "nx": grid[0],
                    "ny": grid[1],
                    "nz": grid[2],
                    "time": float(time),
                    "field": name,
                    "l2": l2,
                    "linf": linf,
                    "relative_l2": l2 / ref_l2 if ref_l2 else None,
                    "relative_linf": linf / ref_linf if ref_linf else None,
                    "native_max_discrepancy": native_discrepancy,
                    "boundary_cell_linf": float(np.max(np.abs(error[[0, -1]]))),
                    "boundary_cell_l2": norms(error[[0, -1]], volume[[0, -1]])[0],
                    "boundary_value_error": boundary_value_error,
                    "boundary_derivative_error": boundary_derivative_error,
                    "max_error_cell": list(map(int, position)),
                }
            )
        psi_phi = fields["phi"] + fields["psi"]
        operator = perpendicular_operator(2 / ds.Bxy**2, psi_phi, ds, mxg=mxg)
        forward = operator.divergence - fields["Vort"]
        kinetic = fields["NVi"] ** 2 / (4 * fields["Ni"]) + fields["NVj"] ** 2 / (
            6 * fields["Nj"]
        )
        partition = fields["Eion"] - 1.5 * fields["Pi"] - kinetic - operator.kperp
        residuals.append(
            {
                "time": float(time),
                "forward_linf": float(abs(forward.isel(selection)).max()),
                "partition_linf": float(abs(partition.isel(selection)).max()),
            }
        )
        require(
            residuals[-1]["forward_linf"] < 5e-10,
            f"Vorticity recovery residual at t={time}: {residuals[-1]['forward_linf']}",
        )
        require(
            residuals[-1]["partition_linf"] < 5e-11,
            f"Energy partition residual at t={time}: {residuals[-1]['partition_linf']}",
        )
    report = {"case": label, "grid": list(grid), "errors": rows, "residuals": residuals}
    (Path(folder) / "errors.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    ds.close()
    return report


def write_scan(folder, reports, *, axis="all", criteria=None):
    """Save all times, orders and full-domain convergence plots; return failures."""
    require(
        not criteria or len(reports) >= 3,
        "Convergence criteria require at least three resolutions; use --no-criteria for smoke runs",
    )
    folder = Path(folder)
    rows = [row.copy() for report in reports for row in report["errors"]]
    for row in rows:
        row.update(l2_order=None, linf_order=None)
    failures = []
    fields = [
        name for name in (*EVOLVED, *RECOVERED) if any(r["field"] == name for r in rows)
    ]
    dimension = {"all": "nx", "x": "nx", "y": "ny", "z": "nz"}[axis]
    for name in fields:
        for time in sorted({row["time"] for row in rows}):
            group = sorted(
                (r for r in rows if r["field"] == name and r["time"] == time),
                key=lambda r: r[dimension],
            )
            for previous, row in pairwise(group):
                for norm in ("l2", "linf"):
                    if row[norm] and previous[norm]:
                        row[norm + "_order"] = float(
                            np.log(previous[norm] / row[norm])
                            / np.log(row[dimension] / previous[dimension])
                        )
            if criteria and len(group) >= 3 and time == max(r["time"] for r in rows):
                for norm in ("l2", "linf"):
                    target = criteria.get(name, criteria.get("default", {})).get(
                        norm + "_min_order"
                    )
                    if target is not None and (
                        group[-1][norm + "_order"] is None
                        or group[-1][norm + "_order"] < target
                    ):
                        failures.append(
                            f"{name}: {norm} order {group[-1][norm + '_order']} < {target}"
                        )
    (folder / "errors.json").write_text(
        json.dumps(rows, indent=2, allow_nan=False) + "\n"
    )
    with (folder / "errors.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(7, 8), constrained_layout=True)
    final = max(r["time"] for r in rows)
    for name in fields:
        group = sorted(
            (r for r in rows if r["field"] == name and r["time"] == final),
            key=lambda r: r[dimension],
        )
        for ax, norm in zip(axes, ("l2", "linf")):
            ax.loglog(
                [float(LENGTH) / r[dimension] for r in group],
                [r[norm] for r in group],
                "o-",
                label=name,
            )
            ax.set(
                xlabel="Mesh spacing (normalised)",
                ylabel=f"{norm} error at t={final:g}",
            )
            ax.grid(True, which="both", alpha=0.3)
    axes[0].legend(ncol=4, fontsize=8)
    for extension in ("pdf", "png"):
        fig.savefig(folder / f"convergence.{extension}", dpi=160)
    plt.close(fig)
    print(
        "case                         field       L2          Linf        L2 order  Linf order"
    )
    for row in rows:
        if row["time"] == final:
            orders = [
                "--" if row[k] is None else f"{row[k]:.3f}"
                for k in ("l2_order", "linf_order")
            ]
            print(
                f"{row['case']:28} {row['field']:7} {row['l2']:11.4e} {row['linf']:11.4e} {orders[0]:>9} {orders[1]:>10}"
            )
    return failures


def compare_cases(first, second):
    """Compare fields directly, avoiding cancellation between error norms."""
    a, ag, grid_a, _ = open_case(first)
    b, bg, grid_b, _ = open_case(second)
    require(
        grid_a == grid_b and ag == bg and np.array_equal(a.t, b.t),
        "Comparison grids/times differ",
    )
    comparisons = {}
    fa, fb = fields_at(a, -1), fields_at(b, -1)
    volume = (
        (a.J * a.dx * a.dy * a.dz)
        .broadcast_like(a.Ni.isel(t=-1))
        .transpose("x", "y", "z")
        .values[ag:-ag]
    )
    for name in fa.keys() & fb.keys():
        error = (fa[name] - fb[name]).transpose("x", "y", "z").values[ag:-ag]
        comparisons[name] = dict(zip(("l2", "linf"), norms(error, volume)))
    a.close()
    b.close()
    return comparisons
