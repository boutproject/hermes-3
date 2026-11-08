"""External powers from resolved options and saved fields, with explicit coverage."""

import numpy as np
import xarray as xr

from .dataset import AuditError
from .fluxes import species_velocity


def coverage(result, missing=()):
    result.attrs.update(
        coverage="incomplete" if missing else "complete",
        missing_terms=sorted(set(missing)),
    )
    return result


def option_bool(options, name, default=False):
    return options.get_bool(name) if name in options else default


def components(run):
    options = run.data.options
    return [
        (
            name.strip(),
            str(options[name.strip()]["type"])
            if name.strip() in options and "type" in options[name.strip()]
            else name.strip(),
        )
        for name in str(options["hermes:components"]).strip(" ()\t").split(",")
    ]


def prescribed_rate(run, section):
    """A scalar or static spatial source in SI, broadcast on the saved state."""
    options, base = run.data.options, xr.ones_like(run.data.Pe)
    key = section + ":source"
    if key not in options:
        if "mesh:file" in options:
            raise AuditError(f"{section} mesh source is not saved")
        return base * 0
    if option_bool(options, section + ":source_time_dependent"):
        raise AuditError(f"{section} time-dependent source prefactor")
    if option_bool(options, section + ":source_only_in_core"):
        raise AuditError(f"{section} core-only source mask")
    try:
        value = options[section].evaluate_scalar("source")
    except (NameError, TypeError):
        if not hasattr(options, "evaluate"):
            raise AuditError(
                f"{section} spatial source coordinates unavailable"
            ) from None
        value = np.asarray(options.evaluate(key))
        shape = tuple(run.data.sizes[d] for d in ("x", "y", "z"))
        try:
            value = np.broadcast_to(value, shape)
        except ValueError:
            raise AuditError(f"{section} source grid differs from saved grid") from None
        value = xr.DataArray(
            value,
            dims=("x", "y", "z"),
            coords={d: run.data[d] for d in ("x", "y", "z")},
        )
    if not bool(np.isfinite(value).all()):
        raise AuditError(f"{section} nonfinite source")
    return base * value


def external_power(run):
    """Per-term external power [W] and particle rates [s^-1].

    Collisional heat exchange is internal and cancels in the plasma reservoir.
    Electron friction needs an evolved electron momentum reservoir; the
    zero-current closure does not supply its mechanical work.
    """
    ds, options = run.data, run.data.options
    fields, missing = {}, []
    volume_scale = run.length_scale**3
    power_scale = run.energy_scale * run.frequency * volume_scale
    zero = run.integrate(xr.zeros_like(ds.Pe))
    fields["power"] = zero.copy()
    if option_bool(options, "solver:mms"):
        missing.append("BOUT++ manufactured forcing")
    for species in run.species:
        name = species.name
        if "evolve_density" in species.components and not option_bool(
            options, "solver:mms"
        ):
            try:
                rate = prescribed_rate(run, "N" + name)
                fields["particles_" + name] = run.integrate(rate) * volume_scale
                if name != "e":
                    v = species_velocity(run, species)
                    normalised = rate / (ds.metadata["Nnorm"] * run.frequency)
                    term = (
                        run.integrate(-0.5 * species.mass * v**2 * normalised)
                        * power_scale
                    )
                    fields["density_kinetic_" + name] = term
                    fields["power"] = fields["power"] + term
            except AuditError as error:
                missing.append(str(error))
        else:
            fields["particles_" + name] = zero.copy()
        if "evolve_pressure" in species.components and not option_bool(
            options, "solver:mms"
        ):
            try:
                term = (
                    1.5 * run.integrate(prescribed_rate(run, "P" + name)) * volume_scale
                )
                fields["pressure_" + name] = term
                fields["power"] = fields["power"] + term
            except AuditError as error:
                missing.append(str(error))
    # Quasineutral electron particle content follows the ion charge sources.
    electron = next(s for s in run.species if s.name == "e")
    if "quasineutral" in electron.components:
        fields["particles_e"] = sum(
            s.charge * fields.get("particles_" + s.name, zero)
            for s in run.species
            if s.name != "e"
        )
    allowed = {
        "evolve_density",
        "evolve_momentum",
        "evolve_pressure",
        "quasineutral",
        "zero_current",
        "braginskii_collisions",
        "braginskii_heat_exchange",
        "sound_speed",
        "evolve_vorticity_ion_energy",
    }
    for name, kinds in components(run):
        for kind in map(str.strip, kinds.split(",")):
            if kind not in allowed:
                missing.append(f"{name}: {kind} power")
        if name in options:
            for flag in (
                "thermal_conduction",
                "low_n_diffuse",
                "low_n_diffuse_perp",
                "low_p_diffuse_perp",
                "low_T_diffuse_perp",
                "damp_p_nt",
            ):
                if option_bool(options[name], flag):
                    missing.append(f"{name}:{flag}")
    result = xr.Dataset(fields).assign_coords(t=ds.t / run.frequency)
    for name in result:
        result[name].attrs.update(
            units="s^-1" if name.startswith("particles_") else "W",
            provenance="resolved physical source options and saved state",
        )
    return coverage(result, missing)
