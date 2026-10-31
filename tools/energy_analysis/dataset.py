"""Loading, species, normalisation and physical-domain contracts."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

from .options import physics_options_file


class AuditError(ValueError):
    """An audit prerequisite is missing or unsupported."""


@dataclass(frozen=True)
class Species:
    name: str
    mass: float
    charge: float
    components: tuple[str, ...] = ()
    density_floor: float = 1e-7


@dataclass
class EnergyRun:
    data: xr.Dataset
    species: tuple[Species, ...]
    model: str
    mxg: int
    myg: int
    energy_scale: float
    length_scale: float
    frequency: float
    source: str = "dataset"
    formulation: str = "total_energy"

    def interior(self, value):
        selection = {}
        for dim, guards in (("x", self.mxg), ("y", self.myg)):
            if dim in value.dims and guards:
                selection[dim] = slice(guards, -guards)
        return value.isel(selection)

    def integrate(self, value):
        volume = self.data.J * self.data.dx * self.data.dy * self.data.dz
        weighted = self.interior(value * volume)
        # Reject NaNs rather than allowing xarray's default sum to hide them.
        return weighted.sum(
            [d for d in ("x", "y", "z") if d in weighted.dims], skipna=False
        )

    def scalar_option(self, name, default=None):
        options = self.data.options[self.model]
        if name not in options:
            if default is None:
                raise AuditError(f"Missing model option {name}")
            return default
        return float(options.evaluate_scalar(name))


def _normalised_view(ds):
    """Convert explicitly annotated SI fields, without mutating the input."""
    result = ds.copy(deep=False)
    result.attrs = dict(ds.attrs)
    result.attrs["metadata"] = dict(ds.attrs.get("metadata", {}))
    coordinates = {}
    for name in list(ds.variables):
        attrs = dict(ds[name].attrs)
        unit_type = attrs.get("units_type")
        if unit_type == "SI":
            conversion = attrs.get("conversion")
            if conversion is None or not np.isfinite(conversion) or conversion <= 0:
                raise AuditError(f"Missing/invalid conversion for SI field {name}")
            value = ds[name] / conversion
            attrs["units_type"] = "hermes"
            value.attrs = attrs
            if name in ds.coords:
                coordinates[name] = value
            else:
                result[name] = value
    # Replace dimensional coordinates last, without aligning physical time
    # labels against the new normalised labels while fields are assigned.
    return result.assign_coords(coordinates)


def _manifest(ds, explicit):
    if explicit is not None:
        species = tuple(s if isinstance(s, Species) else Species(**s) for s in explicit)
    else:
        options = ds.attrs.get("options")
        if options is None:
            raise AuditError(
                "Species discovery requires bd.options or an explicit manifest"
            )
        components = str(options["hermes"]["components"]).strip(" ()\t").split(",")
        species = []
        for component in components:
            name = component.strip()
            if name not in options or "AA" not in options[name]:
                continue
            section = options[name]
            charge = (
                float(section.evaluate_scalar("charge")) if "charge" in section else 0.0
            )
            types = tuple(
                t.strip()
                for t in str(section["type"] if "type" in section else name).split(",")  # noqa: SIM401 - BoutOptions has no get()
            )
            floor = (
                float(section.evaluate_scalar("density_floor"))
                if "density_floor" in section
                else 1e-7
            )
            species.append(
                Species(
                    name, float(section.evaluate_scalar("AA")), charge, types, floor
                )
            )
        species = tuple(species)
    if not species or len({s.name for s in species}) != len(species):
        raise AuditError("Empty or duplicate species manifest")
    for s in species:
        if (
            not np.isfinite(s.mass)
            or s.mass <= 0
            or not np.isfinite(s.charge)
            or not np.isfinite(s.density_floor)
            or s.density_floor <= 0
        ):
            raise AuditError(f"Invalid mass/charge for {s.name}")
        if abs(s.charge) > 1e-5 and "N" + s.name not in ds:
            raise AuditError(f"Missing density for charged species {s.name}")
    if not any(s.name == "e" and s.charge < 0 for s in species):
        raise AuditError("Total ion energy auditing requires electron species e")
    if not any(s.name != "e" and abs(s.charge) > 1e-5 for s in species):
        raise AuditError("No ions found in species manifest")
    return species


def open_run(source, *, species=None, chunks=None):
    """Open a slab run or adapt a fully annotated xhermes Dataset.

    Keep all physical guards, compute in normalised coordinates, and integrate
    only physical cells. Tokamak topology and density-floor corrections are
    deliberately rejected until their reconstruction contracts are validated.
    Load resolved BOUT.settings when available, otherwise supplied BOUT.inp.
    """
    if isinstance(source, xr.Dataset):
        ds = source
        path = "dataset"
    else:
        import xhermes

        path = str(Path(source).resolve())
        options_path = Path(path) / "BOUT.settings"
        if not options_path.is_file():
            options_path = Path(path) / "BOUT.inp"
        with physics_options_file(options_path) as physics_path:
            ds = xhermes.open_hermesdataset(
                datapath=str(Path(path) / "BOUT.dmp.*.nc"),
                inputfilepath=str(physics_path),
                unnormalise=False,
                keep_xboundaries=True,
                keep_yboundaries=True,
                chunks={"t": 1, "y": 1} if chunks is None else chunks,
            )
        if ds.attrs.get("options") is not None and hasattr(ds.options, "filename"):
            ds.options.filename = str(options_path)
    ds = _normalised_view(ds)
    options = ds.attrs.get("options")
    if options is None:
        raise AuditError("Model/boundary reconstruction requires bd.options")
    components = [
        c.strip() for c in str(options["hermes:components"]).strip(" ()\t").split(",")
    ]
    models = [
        (c, str(options[c]["type"]) if c in options and "type" in options[c] else c)
        for c in components
    ]
    models = [
        (c, kind)
        for c, kind in models
        if kind in ("vorticity", "evolve_vorticity_ion_energy")
    ]
    if len(models) != 1:
        raise AuditError(
            "Expected one vorticity or evolve_vorticity_ion_energy component"
        )
    formulation = "pressure" if models[0][1] == "vorticity" else "total_energy"
    manifest = _manifest(ds, species)
    pressure_fields = tuple(
        "P" + s.name for s in manifest if s.name != "e" and abs(s.charge) > 1e-5
    )
    required = (
        "t",
        "Pe",
        *(pressure_fields if formulation == "pressure" else ("Eion", "Tion")),
        "phi",
        "Vort",
        "J",
        "dx",
        "dy",
        "dz",
        "g11",
        "g33",
        "g13",
        "Bxy",
    )
    missing = [name for name in required if name not in ds]
    if missing:
        raise AuditError("Missing inventory fields: " + ", ".join(missing))
    if not {"x", "y", "z"}.issubset(ds.dims):
        raise AuditError(
            "Stage-one slab auditing requires unsqueezed x, y, z dimensions"
        )
    for name in required:
        if ds[name].attrs.get("units_type") != "hermes":
            raise AuditError(
                f"Unknown normalisation for {name}; use xhermes annotated fields"
            )
    metadata = ds.attrs.get("metadata", {})
    rectangle = (
        metadata.get("ixseps1") == 0
        and metadata.get("ixseps2") == 0
        and metadata.get("jyseps1_1") == -1
        and metadata.get("jyseps2_2") == metadata.get("ny", 0) - 1
        and not metadata.get("topology")
    )
    if not rectangle and (
        metadata.get("topology")
        or metadata.get("ixseps1", ds.sizes["x"]) < metadata.get("nx", ds.sizes["x"])
    ):
        raise AuditError("Tokamak boundary masks are not implemented")
    if (
        "mesh:paralleltransform:type" in options
        and str(options["mesh:paralleltransform:type"]) != "identity"
    ):
        raise AuditError("Stage-one auditing requires identity parallel transformation")
    mxg = int(metadata.get("MXG", 0))
    myg = int(metadata.get("MYG", 0)) if metadata.get("keep_yboundaries", False) else 0
    if not metadata.get("keep_xboundaries", False) or mxg < 1:
        raise AuditError("Perpendicular inventory requires retained radial guards")
    if ds.sizes["x"] <= 2 * mxg or ds.sizes["y"] <= 2 * myg:
        raise AuditError("Guard widths leave no physical cells")
    for key in ("rho_s0", "Omega_ci", "Nnorm", "Tnorm"):
        if key not in metadata or not np.isfinite(metadata[key]) or metadata[key] <= 0:
            raise AuditError(f"Missing/invalid normalisation {key}")
    from scipy.constants import elementary_charge

    run = EnergyRun(
        ds,
        manifest,
        models[0][0],
        mxg,
        myg,
        elementary_charge * metadata["Nnorm"] * metadata["Tnorm"],
        metadata["rho_s0"],
        metadata["Omega_ci"],
        path,
        formulation,
    )
    for name in required[1:]:
        if not bool(np.isfinite(run.interior(ds[name])).all().compute()):
            raise AuditError(f"Nonfinite physical field {name}")
    for name in ("J", "dx", "dy", "dz", "Bxy"):
        if bool((run.interior(ds[name]) <= 0).any().compute()):
            raise AuditError(f"Nonpositive geometry field {name}")
    for s in manifest:
        if abs(s.charge) <= 1e-5:
            continue
        density = ds["N" + s.name]
        if density.attrs.get("units_type") != "hermes":
            raise AuditError(f"Unknown normalisation for N{s.name}")
        if not bool(np.isfinite(run.interior(density)).all().compute()):
            raise AuditError(f"Nonfinite density for {s.name}")
        if bool((run.interior(density) <= 20 * s.density_floor).any().compute()):
            raise AuditError(
                f"Density near floor for {s.name}; momentum reconstruction is not validated"
            )
    total_ions = sum(
        ds["N" + s.name] for s in manifest if s.name != "e" and abs(s.charge) > 1e-5
    )
    if bool(
        (run.interior(total_ions) <= 20 * run.scalar_option("density_floor", 1e-7))
        .any()
        .compute()
    ):
        raise AuditError("Total ion density near the common-temperature recovery floor")
    if run.scalar_option("average_atomic_mass") <= 0:
        raise AuditError("Boussinesq average_atomic_mass must be positive")
    return run
