"""Radial wall powers and stored midpoint histories from ordinary snapshots."""

import numpy as np
import xarray as xr

from .dataset import AuditError
from .fluxes import (
    fill_wall_guards,
    polarisation_ac_boundary,
    radial_boundary_mode,
    species_velocity,
)
from .inventory import energy_densities
from .sources import components, coverage, option_bool
from .transfers import polarisation_energy_coefficient
from .transport import (
    flux_dataset,
    fv_xz,
    neighbor,
    solenoidal_reynolds_faces,
    weighted_faces,
)


def absolute_wall_flux(run, transport):
    """Sum magnitudes on both X faces before integrating in time."""
    faces = transport.flux_x.isel(x=[run.mxg - 1, -run.mxg - 1])
    if run.myg:
        faces = faces.isel(y=slice(run.myg, -run.myg))
    return abs(faces).sum(("x", "y", "z"), skipna=False)


def wall_budget(run):
    """Return outward energy/particle fluxes and the measured wall rate.

    Unpinned Reynolds faces exchange weighted quantities at each wall.
    Open chi faces require the standard saved ion_chi output;
    this missing term is reported independently of reconstructable viscosity.
    Time is SI; potential and midpoint rates use Hermes potential units.
    """
    ds = run.data
    zero = run.integrate(xr.zeros_like(ds.Pe))
    fields = {
        "outward": zero.copy(),
        "viscous_outward": zero.copy(),
        "polarisation_outward": zero.copy(),
        "vorticity_rate": zero.copy(),
        "viscous_vorticity_absolute": zero.copy(),
        "polarisation_absolute": zero.copy(),
        "solenoidal_absolute": zero.copy(),
    }
    ions = [s for s in run.species if s.name != "e" and abs(s.charge) > 1e-5]
    fields.update({"particles_" + s.name: zero.copy() for s in ions})
    fields.update({"particles_absolute_" + s.name: zero.copy() for s in ions})
    fields.update({"momentum_" + s.name: zero.copy() for s in ions})
    fields.update({"momentum_absolute_" + s.name: zero.copy() for s in ions})
    missing = []
    periodic = option_bool(ds.options, "periodicX")
    for name, kind in components(run):
        if "sheath" in kind:
            missing.append(f"{name} sheath wall flux")
    for s in run.species:
        if (
            s.name in ds.options
            and option_bool(ds.options[s.name], "bndry_flux", True)
            and not periodic
        ):
            missing.append(f"{s.name} open ExB wall flux")
    if not periodic:
        mode = radial_boundary_mode(run)
        pairs = (("inner", run.mxg, run.mxg - 1), ("outer", -run.mxg - 1, -run.mxg))
        fields["phi_boundary"] = xr.concat(
            [
                ((ds.phi.isel(x=i, drop=True) + ds.phi.isel(x=g, drop=True)) / 2)
                .isel(y=slice(run.myg, ds.sizes["y"] - run.myg))
                .mean(("y", "z"))
                .expand_dims(boundary=[side])
                for side, i, g in pairs
            ],
            "boundary",
        )
        # Rates describe intervals, including their midpoint time, avoiding an
        # invented instantaneous derivative at the first/last output record.
        times = np.asarray(ds.t) / run.frequency
        if len(times) > 1:
            fields["phi_wall_rate"] = xr.DataArray(
                np.diff(fields["phi_boundary"].transpose("boundary", "t"), axis=1)
                / np.diff(times)[None, :],
                dims=("boundary", "interval"),
                coords={
                    "boundary": ["inner", "outer"],
                    "interval": (times[1:] + times[:-1]) / 2,
                },
            )
        scale = run.energy_scale * run.length_scale**3 * run.frequency
        try:
            densities = energy_densities(run)
            Phi = ds.phi + densities.pressure_potential
            if mode == "free_slip" and run.scalar_option("viscosity", 0) > 0:
                nu = run.scalar_option("viscosity") / (
                    run.length_scale**2 * run.frequency
                )
                vort = fill_wall_guards(ds.Vort, run, "zero")
                current = fv_xz(nu, vort, ds, mxg=run.mxg, boundary_flux=True)
                energy = flux_dataset(
                    {
                        d: (Phi + neighbor(Phi, d)) / 2 * current["flux_" + d]
                        for d in ("x", "z")
                    },
                    ds,
                )
                fields["viscous_outward"] = run.integrate(energy.divergence) * scale
                fields["vorticity_rate"] = (
                    run.integrate(current.divergence) * run.frequency
                )
                fields["viscous_vorticity_absolute"] = (
                    absolute_wall_flux(run, current) * run.frequency
                )
            c = run.scalar_option("average_atomic_mass") / ds.Bxy**2
            advected = polarisation_energy_coefficient(run).total
            curl = solenoidal_reynolds_faces(c, ds.phi, Phi, ds)
            ix = xr.DataArray(np.arange(ds.sizes["x"]), dims="x", coords={"x": ds.x})
            wall_mask = (ix == run.mxg - 1) | (ix == ds.sizes["x"] - run.mxg - 1)
            wall_curl = curl["x"].where(wall_mask, 0)
            chi_open = polarisation_ac_boundary(run) == "dirichlet"
            chi = ds.get("ion_chi", ds.get("ddt(Phi)"))
            if chi_open and chi is None:
                missing.append("open polarisation wall flux (chi is not saved)")
            else:
                chi_wall = xr.zeros_like(wall_curl)
                if chi_open:
                    chi_wall = fv_xz(
                        c, chi, ds, mxg=run.mxg, boundary_flux=True
                    ).flux_x.where(wall_mask, 0)

                def wall_transport(weight, *, solenoidal=False):
                    current = -wall_curl if solenoidal else -wall_curl - chi_wall
                    return flux_dataset(
                        weighted_faces({"x": current}, weight, mxg=run.mxg), ds
                    )

                def exchange(transport):
                    # Positive is outward; keep each wall and physical Y row.
                    return xr.concat(
                        [
                            (-transport.flux_x.isel(x=run.mxg - 1, drop=True))
                            .sum("z", skipna=False)
                            .expand_dims(boundary=["inner"]),
                            transport.flux_x.isel(x=-run.mxg - 1, drop=True)
                            .sum("z", skipna=False)
                            .expand_dims(boundary=["outer"]),
                        ],
                        "boundary",
                    ).isel(y=slice(run.myg, ds.sizes["y"] - run.myg))

                pol = wall_transport(advected)
                fields["polarisation_wall_exchange"] = exchange(pol) * scale
                fields["polarisation_outward"] = fields[
                    "polarisation_wall_exchange"
                ].sum(("boundary", "y"))
                fields["polarisation_absolute"] = absolute_wall_flux(run, pol) * scale
                solenoidal = wall_transport(advected, solenoidal=True)
                fields["solenoidal_wall_exchange"] = exchange(solenoidal) * scale
                fields["solenoidal_absolute"] = (
                    absolute_wall_flux(run, solenoidal) * scale
                )
                rho = sum(s.mass * ds["N" + s.name] for s in ions).clip(
                    min=run.scalar_option("average_atomic_mass")
                    * run.scalar_option("density_floor", 1e-7)
                )
                pscale = ds.metadata["Nnorm"] * run.length_scale**3 * run.frequency
                mscale = run.energy_scale * run.length_scale**2
                for species in ions:
                    weight = (
                        species.mass * ds["N" + species.name] / (species.charge * rho)
                    )
                    particle = wall_transport(weight)
                    fields["particles_wall_exchange_" + species.name] = (
                        exchange(particle) * pscale
                    )
                    fields["particles_" + species.name] = (
                        exchange(particle).sum(("boundary", "y")) * pscale
                    )
                    fields["particles_absolute_" + species.name] = (
                        absolute_wall_flux(run, particle) * pscale
                    )
                    momentum = wall_transport(
                        weight * species.mass * species_velocity(run, species)
                    )
                    fields["momentum_wall_exchange_" + species.name] = (
                        exchange(momentum) * mscale
                    )
                    fields["momentum_" + species.name] = (
                        exchange(momentum).sum(("boundary", "y")) * mscale
                    )
                    fields["momentum_absolute_" + species.name] = (
                        absolute_wall_flux(run, momentum) * mscale
                    )
        except AuditError as error:
            missing.append(str(error))
        fields["outward"] = fields["viscous_outward"] + fields["polarisation_outward"]
    fields["particles_e"] = sum(s.charge * fields["particles_" + s.name] for s in ions)
    result = xr.Dataset(fields).assign_coords(t=ds.t / run.frequency)
    for name in result:
        result[name].attrs["units"] = (
            "s^-1"
            if name.startswith("particles_")
            else "N"
            if name.startswith("momentum_")
            else "hermes potential"
            if name == "phi_boundary"
            else "hermes potential/s"
            if name == "phi_wall_rate"
            else "normalised volume/s"
            if name
            in (
                "vorticity_rate",
                "viscous_vorticity_absolute",
            )
            else "W"
        )
    return coverage(result, missing)


def rotation_decay(run):
    """Fit rigid wall rotation and compare with the finite-volume analytic law."""
    from scipy.optimize import curve_fit

    ds = run.data
    history = wall_budget(run).phi_boundary
    difference = np.asarray(
        history.sel(boundary="outer") - history.sel(boundary="inner")
    )
    nx = ds.sizes["x"] - 2 * run.mxg
    length = float((nx * ds.dx / np.sqrt(ds.g11)).values.flat[0]) * run.length_scale
    tau = length / (2 * run.scalar_option("phi_boundary_velocity"))
    fitted, _ = curve_fit(
        lambda time, amplitude, decay: amplitude * (1 - np.exp(-time / decay)),
        np.asarray(history.t),
        difference,
        p0=(difference[-1], tau),
    )
    c = run.scalar_option("average_atomic_mass") / float(ds.Bxy.values.flat[0]) ** 2
    g11 = float(ds.g11.values.flat[0])
    dx = float(ds.dx.values.flat[0])
    omega = ds.Vort.isel(t=0, y=run.myg).mean("z")
    omega1 = float(omega.isel(x=run.mxg + 1) - omega.isel(x=run.mxg)) / dx
    expected = -omega1 * (nx * dx) ** 3 / (12 * c * g11)
    return {
        "tau": tau,
        "fitted_tau": float(fitted[1]),
        "tau_relative_error": abs(float(fitted[1]) / tau - 1),
        "asymptotic_wall_difference": float(fitted[0]),
        "expected_wall_difference": expected,
        "amplitude_error": abs(float(fitted[0]) - expected),
    }


def uniform_drift(run):
    ds = run.data
    budget = wall_budget(run)
    nx = ds.sizes["x"] - 2 * run.mxg
    scalar = lambda name: float(ds[name].values.flat[0])
    c = run.scalar_option("average_atomic_mass") / scalar("Bxy") ** 2
    omega = float(run.interior(ds.Vort).isel(t=0).mean())
    # wall_budget reports potential per SI second, so kappa is divided only
    # by rho_s0 here, rather than by rho_s0 * Omega_ci.
    expected = (
        -run.scalar_option("phi_boundary_velocity")
        / run.length_scale
        * (np.sqrt(scalar("g11")) * omega * nx * scalar("dx") / (2 * c * scalar("g11")))
    )
    return float(abs(budget.phi_wall_rate / expected - 1).max())
