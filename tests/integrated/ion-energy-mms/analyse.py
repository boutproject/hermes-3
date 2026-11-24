"""Continuum errors and convergence plots, independent of simulation launching."""

import manufactured as exact
import numpy as np
import pandas as pd
import xarray as xr


def fields(ds):
    result = {name: ds[name] for name in exact.EVOLVED}
    result.update(
        phi=ds.phi,
        psi=ds.ion_pressure_potential,
        chi=ds.ion_chi,
        Tion=ds.Tion,
        pi=ds.Ni * ds.Tion,
        pj=ds.Nj * ds.Tion,
        Pi=(ds.Ni + ds.Nj) * ds.Tion,
        Te=ds.Te,
    )
    return result


def norms(ds, error):
    error = error.assign_attrs(cell_location="CELL_CENTRE", direction_y="Standard")
    volume = ds.hermes.integrate_midpoints(
        xr.ones_like(error).assign_attrs(error.attrs)
    )
    return {
        "l2": float(np.sqrt(ds.hermes.integrate_midpoints(error**2) / volume)),
        "linf": float(abs(error).max()),
    }


def errors(ds):
    """Final-time global norms and norms over the fixed central radial half."""
    time = float(ds.t.isel(t=-1))
    ds = ds.isel(t=-1, drop=True)
    length = float(exact.LENGTH)
    x = length * (np.arange(ds.sizes["x"]) + 0.5) / ds.sizes["x"]
    y = length * (np.arange(ds.sizes["y"]) + 0.5) / ds.sizes["y"]
    z = length * np.arange(ds.sizes["z"]) / ds.sizes["z"]
    interior = np.flatnonzero((x >= length / 4) & (x <= 3 * length / 4))
    result = {}
    for name, field in fields(ds).items():
        reference = exact.evaluate(
            name, x[:, None, None], y[None, :, None], z[None, None, :], time
        )
        error = field - xr.DataArray(
            reference,
            dims=("x", "y", "z"),
            coords={dim: ds[dim] for dim in ("x", "y", "z")},
        )
        result[name] = norms(ds, error)
        result[name].update(
            {
                key + "_interior": value
                for key, value in norms(
                    ds.isel(x=interior), error.isel(x=interior)
                ).items()
            }
        )
    return result


def check_recovery(ds):
    recovered = 1.5 * (ds.Ni + ds.Nj) * ds.Tion + ds.ion_KE_perp + ds.ion_KE_parallel
    assert float(abs(recovered - ds.Eion).max()) < 1e-10
    assert float(abs(ds.ion_auxiliary_residual).max()) < 1e-10
    for name in (
        "parallel_current",
        "DivJpar",
        "DivJReynolds",
        "DivJviscous",
        "DivJpol",
    ):
        assert float(abs(ds[name]).max()) > 1e-6, name


def save_convergence(output, reports):
    import matplotlib.pyplot as plt

    rows = [
        dict(nx=size, field=name, **values)
        for size, report in reports.items()
        for name, values in report.items()
    ]
    pd.DataFrame(rows).to_csv(output / "errors.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for name in next(iter(reports.values())):
        for axis, norm in zip(axes, ("l2_interior", "linf_interior")):
            axis.loglog(
                list(reports),
                [report[name][norm] for report in reports.values()],
                "o-",
                label=name,
            )
            axis.set(xlabel="Radial cells", ylabel=norm.replace("_", " "))
    axes[-1].legend(fontsize="small", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.savefig(output / "convergence.pdf")
