"""Limited snapshot reconstruction of internal energy transfer in slab cases."""

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
