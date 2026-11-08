"""Limited snapshot reconstruction of internal energy transfer in slab cases."""

import numpy as np
import xarray as xr

from .dataset import AuditError
from .inventory import _kinetic
from .operators import spectral_ddz


def electron_curvature_work(run):
    """Electron thermal ExB work, normalised; supports analytic Z-only curvature.

    Requires an identity slab, scalar input curvature and an unfiltered FFT
    DDZ. Grid-file curvature and other derivative schemes need their own
    reconstruction. No optional diagnostic is used.
    """
    options = run.data.options
    if "mesh:file" in options:
        raise AuditError("Grid-file curvature reconstruction is unavailable")
    for component in ("x", "y", "z"):
        key = "bxcv" + component
        if key not in options["mesh"]:
            raise AuditError(f"Explicit scalar mesh:{key} is required")
    if any(options["mesh"].evaluate_scalar("bxcv" + c) != 0 for c in ("x", "y")):
        raise AuditError("Curvature work currently supports Z-only curvature")
    if (
        "mesh:ddz:first" not in options
        or str(options["mesh:ddz:first"]).upper() != "FFT"
    ):
        raise AuditError("Curvature work requires explicit FFT DDZ")
    if (
        "mesh:ddz:fft_filter" not in options
        or options["mesh:ddz"].evaluate_scalar("fft_filter") != 0
    ):
        raise AuditError("Curvature work requires unfiltered FFT DDZ")
    curl_z = (
        2
        * run.length_scale**2
        * options["mesh"].evaluate_scalar("bxcvz")
        / run.data.Bxy
    )
    work = -run.data.Pe * curl_z * spectral_ddz(run.data.phi, run.data.dz)
    work.attrs = {
        "units_type": "hermes",
        "provenance": "reconstructed Z curvature and spectral DDZ from snapshots",
    }
    return work


def integrated_power(run, density):
    """Convert a normalised energy RHS density to a scalar power history in W."""
    if density.attrs.get("units_type") != "hermes":
        raise AuditError("Power reconstruction requires a normalised RHS density")
    result = run.integrate(density) * (
        run.energy_scale * run.length_scale**3 * run.frequency
    )
    result = result.assign_coords(t=run.data.t / run.frequency)
    result.t.attrs = {"units": "s"}
    result.attrs = {
        "units": "W",
        "provenance": density.attrs.get("provenance", "supplied"),
    }
    return result


def collisional_heat_exchange(run):
    """Reconstruct Braginskii electron-ion heat transfer in a common-T slab.

    Ion-ion heat transfer is exactly zero with a common ion temperature.
    Frequencies use the component's NRL Coulomb logarithm and SI constants;
    the returned per-species powers cancel in the plasma reservoir.
    """
    from scipy.constants import electron_mass as me
    from scipy.constants import elementary_charge as qe
    from scipy.constants import epsilon_0
    from scipy.constants import proton_mass as mp

    from .inventory import soft_floor
    from .sources import components, coverage, option_bool

    ds = run.data
    zero = xr.zeros_like(ds.Pe)
    fields = {s.name: zero.copy() for s in run.species}
    collision = next(
        (name for name, kind in components(run) if kind == "braginskii_collisions"),
        None,
    )
    exchange = any(kind == "braginskii_heat_exchange" for _, kind in components(run))
    if not collision or not exchange:
        return coverage(xr.Dataset())
    electron = next(s for s in run.species if s.name == "e")
    if option_bool(ds.options[collision], "electron_ion", True):
        temperature = ds.metadata["Tnorm"]
        density = ds.metadata["Nnorm"]
        te = ds.Pe / ds.Ne * temperature
        ti = ds.Tion * temperature
        ne = ds.Ne * density
        multiplier = (
            float(ds.options[collision].evaluate_scalar("ei_multiplier"))
            if "ei_multiplier" in ds.options[collision]
            else 1.0
        )
        for ion in run.species:
            if ion.name == "e" or ion.charge <= 0:
                continue
            ni = ds["N" + ion.name] * density
            ratio = me / (mp * ion.mass)
            log = xr.where(
                (te < 0.1) | (ni < 1e10) | (ne < 1e10),
                10,
                xr.where(
                    te < ti * ratio,
                    23
                    - 0.5 * np.log(ni)
                    + 1.5 * np.log(ti)
                    - np.log(ion.charge**2 * ion.mass),
                    xr.where(
                        te < np.exp(2) * ion.charge**2,
                        30 - 0.5 * np.log(ne) - np.log(ion.charge) + 1.5 * np.log(te),
                        31 - 0.5 * np.log(ne) + np.log(te),
                    ),
                ),
            )
            ve2 = 2 * soft_floor(te, 0.1) * qe / me
            vi2 = 2 * soft_floor(ti, 0.1) * qe / (mp * ion.mass)
            frequency = (
                qe**4
                * ion.charge**2
                * ni.clip(min=0)
                * soft_floor(log, 1)
                * (1 + ratio)
                / (3 * (np.pi * (ve2 + vi2)) ** 1.5 * (epsilon_0 * me) ** 2)
                * multiplier
                / run.frequency
            )
            transfer = (
                3
                * electron.mass
                / (electron.mass + ion.mass)
                * frequency
                * ds.Ne
                * (ds.Tion - ds.Pe / ds.Ne)
            )
            fields["e"] = fields["e"] + transfer
            fields[ion.name] = fields[ion.name] - transfer
    result = xr.Dataset()
    for name, field in fields.items():
        field.attrs.update(
            units_type="hermes", provenance="reconstructed Braginskii heat exchange"
        )
        result[name] = integrated_power(run, field)
    result["net"] = sum(result.values())
    result["net"].attrs["units"] = "W"
    missing = (
        ["neutral collisional heat exchange"]
        if any(s.charge == 0 for s in run.species)
        else []
    )
    return coverage(result, missing)


def polarisation_energy_coefficient(run):
    """Multispecies E_adv, including parallel kinetic energy and mass floor.

    This reconstructs the transported coefficient, not the polarisation flux.
    """
    ions = [s for s in run.species if s.name != "e" and abs(s.charge) > 1e-5]
    rho = sum(s.mass * run.data["N" + s.name] for s in ions)
    floor = run.scalar_option("density_floor", 1e-7) * run.scalar_option(
        "average_atomic_mass"
    )
    rho = rho.clip(min=floor)
    enthalpy = sum(
        s.mass / s.charge * 2.5 * run.data["N" + s.name] * run.data.Tion / rho
        for s in ions
    )
    kinetic = sum(s.mass / s.charge * _kinetic(run, s) / rho for s in ions)
    result = xr.Dataset(
        {
            "enthalpy": enthalpy,
            "parallel_kinetic": kinetic,
            "total": run.data.phi + enthalpy + kinetic,
        }
    )
    result.attrs = {
        "units_type": "hermes",
        "provenance": "reconstructed multispecies energy coefficient",
    }
    return result
