"""Independent flux/source audits for the implemented common-temperature slab model.

The supported snapshot contract is deliberately narrower than the inventory
contract: uniform orthogonal metrics, identity transform, unfiltered FFT Z,
C2 parallel derivatives, MC parallel limiter, and no auxiliary dissipation.
No diagnostic is used to reconstruct an RHS or solve for the polarisation flow.
"""

from dataclasses import dataclass

import numpy as np
import xarray as xr

from .dataset import AuditError
from .inventory import _kinetic, energy_densities, soft_floor
from .operators import spectral_ddz
from .transfers import (
    electron_curvature_work,
    polarisation_energy_coefficient,
)
from .transport import (
    arakawa,
    bout_perp,
    central_parallel,
    flux_dataset,
    fv_xz,
    neighbor,
    numerical_heating_defect,
    numerical_parallel_heating,
    parallel_flux,
    parallel_kinetic_flux,
    spectral_d2z,
    xppm,
)


@dataclass(frozen=True)
class PrescribedSources:
    """External SI rates per species, with explicit coverage and provenance.

    Each species must have an entry (even for zero); rates may be scalars or
    labelled fields on the run's normalised coordinates. Thermal is W/m³,
    density is m⁻³/s, momentum is N/m³. Internal polarisation sources are excluded.
    Electron thermal heating is supplied under 'e'.
    """

    thermal: dict
    density: dict
    momentum: dict
    provenance: str


def species_velocity(run, species):
    ds = run.data
    if "NV" + species.name in ds:
        return ds["NV" + species.name] / (
            species.mass * soft_floor(ds["N" + species.name], species.density_floor)
        )
    if "V" + species.name in ds:
        return xr.ones_like(ds["N" + species.name]) * ds["V" + species.name]
    if set(species.components) & {
        "evolve_momentum",
        "zero_current",
        "fixed_velocity",
        "electromagnetic",
    }:
        raise AuditError(
            f"Missing saved momentum/velocity for flowing species {species.name}"
        )
    return xr.zeros_like(ds["N" + species.name])


def reconstruct_sources(run, sources):
    """Thermal plus momentum/particle kinetic conversion from prescribed SI rates."""
    if not sources.provenance:
        raise AuditError("Prescribed sources require provenance")
    names = {s.name for s in run.species}
    if any(
        set(rates) != names
        for rates in (sources.thermal, sources.density, sources.momentum)
    ):
        raise AuditError("Explicit source coverage is required for every species")
    from scipy.constants import proton_mass

    nscale = run.data.metadata["Nnorm"] * run.frequency
    mscale = (
        proton_mass * run.data.metadata["Nnorm"] * run.length_scale * run.frequency**2
    )
    qscale = run.energy_scale * run.frequency
    fields = {}
    ion = xr.zeros_like(run.data.Eion)
    for s in run.species:
        n = xr.ones_like(run.data.Eion) * sources.density[s.name] / nscale
        momentum_scale = (
            run.data["NV" + s.name].attrs["conversion"] * run.frequency
            if "Cs0" in run.data.metadata and "NV" + s.name in run.data
            else mscale
        )
        m = xr.ones_like(run.data.Eion) * sources.momentum[s.name] / momentum_scale
        q = xr.ones_like(run.data.Eion) * sources.thermal[s.name] / qscale
        v = species_velocity(run, s)
        total = q + v * m - 0.5 * s.mass * v**2 * n
        fields.update(
            {
                "density_" + s.name: n,
                "momentum_" + s.name: m,
                "thermal_" + s.name: q,
                "total_" + s.name: total,
            }
        )
        if s.name != "e":
            ion += total
    fields["Eion_src"] = ion
    fields["plasma_src"] = ion + fields["total_e"]
    result = xr.Dataset(fields)
    for f in result:
        result[f].attrs = {"units_type": "hermes", "provenance": sources.provenance}
    return result


def slab_contract(run):
    ds, options, meta = run.data, run.data.options, run.data.metadata
    for name in ("J", "dx", "dy", "dz", "Bxy", "g11", "g33", "g_22"):
        if name not in ds:
            raise AuditError(f"Flux reconstruction requires saved {name}")
        values = np.asarray(ds[name])
        if not np.isfinite(values).all() or not np.allclose(
            values, values.flat[0], rtol=1e-12, atol=1e-14
        ):
            raise AuditError(f"Flux reconstruction requires uniform saved {name}")
    for name in ("g13", "g23", "g_23"):
        if name not in ds or float(abs(ds[name]).max()) > 1e-14:
            raise AuditError(f"Flux reconstruction requires orthogonal saved {name}=0")
    if str(options["mesh:paralleltransform:type"]) != "identity":
        raise AuditError("Flux reconstruction requires identity parallel transform")
    if (
        str(options["mesh:ddz:first"]).upper() != "FFT"
        or options["mesh:ddz"].evaluate_scalar("fft_filter") != 0
    ):
        raise AuditError("Flux reconstruction requires unfiltered FFT DDZ")
    for key in ("mesh:ddy:first", "ddy:first"):
        if key in options and str(options[key]).upper() != "C2":
            raise AuditError("Parallel work reconstruction requires C2 DDY")
    if str(meta.get("HERMES_SLOPE_LIMITER", "")) != "MC":
        raise AuditError("Parallel flux reconstruction requires saved MC limiter")
    for key in ("viscosity", "hyper_z"):
        if run.scalar_option(key, 0) > 0:
            raise AuditError(f"Unaccounted ion-energy {key}")
    if options[run.model].get_bool("phi_boundary_relax") or options[run.model].get_bool(
        "poloidal_flows"
    ):
        raise AuditError(
            "Potential relaxation/poloidal flow reconstruction is unavailable"
        )
    allowed = {
        "evolve_density",
        "evolve_pressure",
        "evolve_momentum",
        "quasineutral",
        "zero_current",
        "fixed_velocity",
        "test_energy_sources",
    }
    for s in run.species:
        if abs(s.charge) <= 1e-5:
            raise AuditError(
                "Neutral-species energy transport is outside flux coverage"
            )
        if set(s.components) - allowed:
            raise AuditError(f"Unaccounted physics components for {s.name}")
        section = options[s.name]
        for key in (
            "thermal_conduction",
            "poloidal_flows",
            "low_n_diffuse",
            "low_n_diffuse_perp",
            "low_p_diffuse_perp",
            "low_T_diffuse_perp",
            "damp_p_nt",
            "evolve_log",
            "p_div_v",
        ):
            if key in section and section.get_bool(key):
                raise AuditError(f"Unaccounted {s.name}:{key}")
        for key in ("hyper_z", "hyper_z_T"):
            if key in section and section.evaluate_scalar(key) > 0:
                raise AuditError(f"Unaccounted {s.name}:{key}")
        if "evolve_momentum" in s.components and not section.get_bool(
            "use_div_par_fvv"
        ):
            raise AuditError("Parallel momentum reconstruction requires Div_par_fvv")
    components = {
        c.strip() for c in str(options["hermes:components"]).strip(" ()\t").split(",")
    }
    if components != {run.model, *(s.name for s in run.species)}:
        raise AuditError("Additional components leave source/flux coverage incomplete")
    # Exactly a core-periodic or all-open rectangle, without branch connections.
    periodic_y = meta.get("ixseps1") == meta.get("nx")
    rectangle = meta.get("ixseps1") == 0 and meta.get("ixseps2") == 0
    if (
        not (periodic_y or rectangle)
        or meta.get("jyseps1_1") != -1
        or meta.get("jyseps2_2") != meta.get("ny") - 1
    ):
        raise AuditError(
            "Flux reconstruction requires simple periodic/open rectangular Y topology"
        )
    if not periodic_y and run.myg < 2:
        raise AuditError("Open parallel fluxes require two retained Y halos")
    return options.get_bool("periodicX"), periodic_y


def fill_transport_guards(field, run, periodic_x, periodic_y):
    """Communication or homogeneous neumann_o2 for internally computed fields.

    BOUT's second Neumann guard mirrors the second interior cell.
    """
    result = field.copy()
    for dim, width, periodic in (
        ("x", run.mxg, periodic_x),
        ("y", run.myg, periodic_y),
    ):
        if not width:
            continue
        if width != 2:
            raise AuditError(
                "Internal transport guard reconstruction requires width two"
            )
        n = result.sizes[dim]
        for i in range(width):
            low = n - 2 * width + i if periodic else 2 * width - 1 - i
            high = width + i if periodic else n - width - 1 - i
            result[{dim: i}] = result.isel({dim: low})
            result[{dim: n - width + i}] = result.isel({dim: high})
    return result


def solve_polarisation(run, rhs, periodic_x, periodic_y):
    """Independent sparse FV solve, gauge-fixed for periodic X.

    Nonperiodic X uses homogeneous Dirichlet at half-cell faces (default
    Laplace flags). The returned transport field has Neumann radial guards,
    as in the component; its post-boundary forward residual can differ.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import splu

    ds = run.data
    # Nondefault inversion boundary flags or coefficients need explicit support.
    section = ds.options[run.model]
    for key in (
        "laplacian:inner_boundary_flags",
        "laplacian:outer_boundary_flags",
        "laplacian:global_flags",
    ):
        if key in section and section.evaluate_scalar(key) != 0:
            raise AuditError(
                "Polarisation solve supports only default Laplace boundary flags"
            )
    nx, nz = ds.sizes["x"] - 2 * run.mxg, ds.sizes["z"]
    c = run.scalar_option("average_atomic_mass") / float(ds.Bxy.values.flat[0]) ** 2
    ax = c * float(ds.g11.values.flat[0]) / float(ds.dx.values.flat[0]) ** 2
    az = c * float(ds.g33.values.flat[0]) / float(ds.dz.values.flat[0]) ** 2
    rows, cols, entries = [], [], []
    for x in range(nx):
        for z in range(nz):
            i = x * nz + z
            diagonal = -2 * ax - 2 * az
            for xx in (x - 1, x + 1):
                if periodic_x or 0 <= xx < nx:
                    rows.append(i)
                    cols.append((xx % nx) * nz + z)
                    entries.append(ax)
                else:
                    diagonal -= ax  # reflected negative guard at the face
            for zz in (z - 1, z + 1):
                rows.append(i)
                cols.append(x * nz + (zz % nz))
                entries.append(az)
            rows.append(i)
            cols.append(i)
            entries.append(diagonal)
    matrix = coo_matrix((entries, (rows, cols)), shape=(nx * nz, nx * nz)).tocsc()
    if periodic_x:
        matrix = matrix.tolil()
        matrix[0, :] = 0
        matrix[0, 0] = 1
        matrix = matrix.tocsc()
    factor = splu(matrix)

    def solve(values):
        b = np.asarray(values).ravel().copy()
        if not np.isfinite(b).all():
            raise AuditError("Polarisation RHS contains nonfinite physical values")
        if periodic_x:
            if abs(b.mean()) > 1e-11 * max(1.0, abs(b).max()):
                raise AuditError(
                    "Periodic polarisation RHS has incompatible nonzero mean"
                )
            b[0] = 0
        solved = factor.solve(b).reshape(nx, nz)
        return solved - solved.mean() if periodic_x else solved

    interior = run.interior(rhs).reset_coords(drop=True)
    if interior.chunks is not None:
        interior = interior.chunk({"x": -1, "z": -1})
    solved = xr.apply_ufunc(
        solve,
        interior,
        input_core_dims=[["x", "z"]],
        output_core_dims=[["x", "z"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float],
    )
    # Reindex introduces only undefined guard placeholders; communication or
    # the known Neumann rule fills them lazily. No full-time 4D allocation.
    result = solved.reindex(x=rhs.x, y=rhs.y).transpose(*rhs.dims)
    return fill_transport_guards(result, run, periodic_x, periodic_y)


def reconstruct_fluxes(run, *, sources=None):
    """Return supported energy/current divergences, face fluxes and species RHS.

    Without prescribed sources fluxes remain usable, but external sources and
    full RHS identities are explicitly unavailable. Unknown schemes raise
    AuditError; audit_report converts that to incomplete coverage.
    """
    periodic_x, periodic_y = slab_contract(run)
    ds = run.data
    densities = energy_densities(run)
    psi, phi = densities.pressure_potential, ds.phi
    Phi = phi + psi
    c = run.scalar_option("average_atomic_mass") / ds.Bxy**2
    pi = sum(ds["N" + s.name] * ds.Tion for s in run.species if s.name != "e")
    bracket = fill_transport_guards(
        arakawa(phi, psi, ds.dx, ds.dz), run, periodic_x, periodic_y
    )
    lap = ds.g11 * (
        neighbor(phi, "x") - 2 * phi + neighbor(phi, "x", -1)
    ) / ds.dx**2 + ds.g33 * spectral_d2z(phi, ds.dz)
    delp = fill_transport_guards(0.5 * c * lap, run, periodic_x, periodic_y)
    families = {}

    def combine(*terms):
        fluxes = {
            d: sum(
                sign * term["flux_" + d] for sign, term in terms if "flux_" + d in term
            )
            for d in ("x", "y", "z")
            if any("flux_" + d in term for _, term in terms)
        }
        return flux_dataset(fluxes, ds)

    reynolds = combine(
        (-1, xppm(0.5 * ds.Vort, phi, ds)),
        (-1, bout_perp(0.5 * c, bracket, ds, myg=run.myg)),
        (-1, xppm(delp, Phi, ds)),
    )
    curl = (
        2 * run.length_scale**2 * ds.options["mesh"].evaluate_scalar("bxcvz") / ds.Bxy
    )
    # electron_curvature_work validates scalar, Z-only curvature and DDZ settings.
    work = electron_curvature_work(run)
    dia = spectral_ddz((ds.Pe + pi) * curl, ds.dz)
    dphi = solve_polarisation(run, reynolds.divergence + dia, periodic_x, periodic_y)

    def pol(weight):
        return combine(
            (
                -1,
                fv_xz(weight * c, dphi, ds, boundary_flux=False, periodic_x=periodic_x),
            ),
            (1, xppm(0.5 * ds.Vort, -phi, ds, weight=weight)),
            (
                -1,
                fv_xz(
                    0.5 * weight * c,
                    bracket,
                    ds,
                    boundary_flux=False,
                    periodic_x=periodic_x,
                ),
            ),
            (1, xppm(delp, -Phi, ds, weight=weight)),
        )

    families["DivEExB"] = xppm(1.5 * pi + densities.ion_parallel, phi, ds)
    families["DivEReynolds"] = bout_perp(0.5 * c * bracket, Phi, ds, myg=run.myg)
    families["DivEpol"] = pol(polarisation_energy_coefficient(run).total)
    families["DivJpol"] = pol(xr.ones_like(phi))
    families["DivJReynolds"] = reynolds
    result = {
        "ddt(Phi)": dphi,
        "ddt(Vort)": reynolds.divergence + dia,
        "electron_exb_work": work,
        "arakawa": bracket,
    }
    jion = xr.zeros_like(phi)
    je = xr.zeros_like(phi)
    parallel = []
    rho = sum(s.mass * ds["N" + s.name] for s in run.species if s.name != "e").clip(
        min=run.scalar_option("average_atomic_mass")
        * run.scalar_option("density_floor", 1e-7)
    )
    source_fields = reconstruct_sources(run, sources) if sources is not None else None
    for s in run.species:
        name, n = s.name, ds["N" + s.name]
        v = species_velocity(run, s)
        p = ds.Pe if name == "e" else n * ds.Tion
        wave = np.sqrt(p / (n * s.mass))
        if name == "e":
            je = s.charge * n * v
        else:
            jion += s.charge * n * v
            parallel.append(
                (
                    1,
                    parallel_flux(
                        2.5 * p + _kinetic(run, s),
                        v,
                        wave,
                        ds,
                        myg=run.myg,
                        periodic_y=periodic_y,
                    ),
                )
            )
            families["particle_pol_" + name] = pol(s.mass * n / (s.charge * rho))
            families["momentum_pol_" + name] = pol(
                s.mass * s.mass * n * v / (s.charge * rho)
            )
        if source_fields is not None:
            density_rhs = (
                -xppm(
                    n,
                    phi,
                    ds,
                    boundary_flux=ds.options[name].get_bool("bndry_flux")
                    if "bndry_flux" in ds.options[name]
                    else False,
                ).divergence
                - parallel_flux(
                    n, v, wave, ds, myg=run.myg, periodic_y=periodic_y, fixflux=False
                ).divergence
                + source_fields["density_" + name]
            )
            if name != "e":
                density_rhs -= families["particle_pol_" + name].divergence
            result["ddt(N" + name + ")"] = density_rhs
            if "evolve_momentum" in s.components:
                rhs = (
                    -xppm(
                        s.mass * n * v,
                        phi,
                        ds,
                        boundary_flux=ds.options[name].get_bool("bndry_flux")
                        if "bndry_flux" in ds.options[name]
                        else False,
                    ).divergence
                    - s.mass
                    * parallel_flux(
                        n,
                        v,
                        wave,
                        ds,
                        myg=run.myg,
                        periodic_y=periodic_y,
                        momentum=True,
                    ).divergence
                    - central_parallel(p, ds)
                    - s.charge * n * central_parallel(phi, ds)
                    + source_fields["momentum_" + name]
                )
                if name != "e":
                    rhs -= families["momentum_pol_" + name].divergence
                result["ddt(NV" + name + ")"] = rhs
    families["DivEparallel"] = combine(*parallel)

    def electric_flux(field_z, field_y):
        divz = spectral_ddz(field_z, ds.dz)
        # Exact face antiderivative of the spectral divergence. The zero flux
        # gauge is arbitrary on periodic Z and does not change regional budgets.
        fz = (divz * ds.J * ds.dx * ds.dy * ds.dz).cumsum("z", skipna=False)
        yp = lambda v: neighbor(v, "y")
        area = ds.J * ds.dx * ds.dz / np.sqrt(ds.g_22)
        fy = (field_y + yp(field_y)) / 2 * area if run.myg else xr.zeros_like(field_y)
        return flux_dataset({"z": fz, "y": fy}, ds)

    families["DivEelectric"] = electric_flux(phi * pi * curl, phi * jion)
    result["phi_Div_Je"] = phi * (
        spectral_ddz(ds.Pe * curl, ds.dz) + central_parallel(je, ds)
    )
    if source_fields is not None:
        result["Eion_src"] = source_fields.Eion_src
        result["ddt(Eion)"] = (
            -sum(
                families[n].divergence
                for n in (
                    "DivEExB",
                    "DivEReynolds",
                    "DivEpol",
                    "DivEparallel",
                    "DivEelectric",
                )
            )
            - result["phi_Div_Je"]
            + source_fields.Eion_src
        )
        e = next(s for s in run.species if s.name == "e")
        ve = species_velocity(run, e)
        we = np.sqrt(ds.Pe / (ds.Ne * e.mass))
        families["electron_ExB"] = xppm(
            1.5 * ds.Pe, phi, ds, boundary_flux=ds.options["e"].get_bool("bndry_flux")
        )
        families["electron_parallel_thermal"] = parallel_flux(
            2.5 * ds.Pe, ve, we, ds, myg=run.myg, periodic_y=periodic_y
        )
        heating = (
            e.mass
            * numerical_parallel_heating(soft_floor(ds.Ne, e.density_floor), ve, we, ds)
            if ds.options["e"].get_bool("numerical_viscous_heating") and run.myg
            else xr.zeros_like(phi)
        )
        result["ddt(Pe)"] = (2 / 3) * (
            -families["electron_ExB"].divergence
            - families["electron_parallel_thermal"].divergence
            + ve * central_parallel(ds.Pe, ds)
            + work
            + source_fields.thermal_e
            + heating
        )
        if "quasineutral" in e.components:
            result["ddt(Ne)"] = sum(
                s.charge * result["ddt(N" + s.name + ")"]
                for s in run.species
                if s.name != "e"
            )
        for q in run.species:
            if "quasineutral" in q.components:
                result["ddt(N" + q.name + ")"] = (
                    -sum(
                        other.charge * result["ddt(N" + other.name + ")"]
                        for other in run.species
                        if other.name != q.name
                    )
                    / q.charge
                )
        result["electron_numerical_heating"] = heating
        result["electron_heating_defect"] = (
            e.mass * numerical_heating_defect(ds.Ne, ve, we, ds)
            if ds.options["e"].get_bool("numerical_viscous_heating") and run.myg
            else xr.zeros_like(phi)
        )
        if "evolve_momentum" in e.components and periodic_y:
            families["electron_parallel_kinetic"] = (
                parallel_kinetic_flux(ds.Ne, ve, we, ds, myg=run.myg) * e.mass
            )
        result["electron_kinetic_rhs"] = (
            ve * result["ddt(NVe)"] - 0.5 * e.mass * ve**2 * result["ddt(Ne)"]
            if "evolve_momentum" in e.components
            else xr.zeros_like(phi)
        )
    for name, family in families.items():
        result[name] = family.divergence
        result.update(
            {name + "_" + f: family[f] for f in family if f.startswith("flux_")}
        )
    output = xr.Dataset(result)
    for name in output:
        output[name].attrs = {"units_type": "hermes", "provenance": "reconstructed"}
    for family in families:
        quantity = (
            "energy"
            if family.startswith(("DivE", "electron_ExB", "electron_parallel"))
            else "current"
            if family.startswith("DivJ")
            else "particles"
            if family.startswith("particle_pol")
            else "momentum"
        )
        for name in output:
            if name == family or name.startswith(family + "_flux_"):
                output[name].attrs["quantity"] = quantity
    output.attrs = {
        "flux_orientation": "increasing coordinate; lower-cell index",
        "flux_measure": "area integrated",
        "source_coverage": sources.provenance if sources else "unavailable",
        "geometry_coverage": "uniform orthogonal identity slab",
        "periodic_x": periodic_x,
        "periodic_y": periodic_y,
        "limitations": "XPPM closes the periodic radial seam when bndry_flux=false; BOUT Reynolds FV omits Z when MYG=0",
    }
    return output


def net_face_flux(run, fluxes, family, *, selection=None):
    """Normalised net outward face flux for an energy/current/species family.

    Selection contains contiguous cell slices in X/Y/Z; omit to select all
    physical cells. Internal and periodic interfaces cancel by construction.
    """
    selection = {} if selection is None else selection
    slices = {
        d: selection.get(d, slice(g, run.data.sizes[d] - g))
        for d, g in (("x", run.mxg), ("y", run.myg), ("z", 0))
    }
    terms = []
    for dim in ("x", "y", "z"):
        name = family + "_flux_" + dim
        if name not in fluxes:
            continue
        face = fluxes[name]
        section = slices[dim]
        if (
            section.step not in (None, 1)
            or section.start is None
            or section.stop is None
            or not 0 <= section.start < section.stop <= face.sizes[dim]
        ):
            raise AuditError(
                "Subvolume budgets require explicit contiguous cell slices"
            )
        periodic_without_guards = (
            dim == "y" and run.myg == 0 and fluxes.attrs.get("periodic_y")
        )
        low = (
            (section.start - 1) % face.sizes[dim]
            if dim == "z" or periodic_without_guards
            else section.start - 1
        )
        if low < 0:
            raise AuditError("Flux budget requires the lower saved face halo")
        diff = face.isel({dim: section.stop - 1}) - face.isel({dim: low})
        diff = diff.isel({d: s for d, s in slices.items() if d in diff.dims})
        terms.append(
            diff.sum([d for d in ("x", "y", "z") if d in diff.dims], skipna=False)
        )
    if not terms:
        raise AuditError(f"No reconstructed faces for {family}")
    result = sum(terms).reset_coords(drop=True)
    result.attrs = {
        "units_type": "hermes",
        "provenance": "reconstructed area-integrated outward face flux",
        "quantity": fluxes[family].attrs.get("quantity", "unspecified"),
    }
    return result


def flux_power(run, fluxes, family, *, selection=None):
    """Net outward energy power in W; current/species fluxes are rejected.

    Use net_face_flux for normalised currents or particle/momentum transport.
    Energy faces are already area-integrated; no additional area is applied.
    """
    if family not in fluxes or fluxes[family].attrs.get("quantity") != "energy":
        raise AuditError(
            "flux_power requires an energy-flux family; use net_face_flux for currents or species fluxes"
        )
    result = (
        net_face_flux(run, fluxes, family, selection=selection)
        * run.energy_scale
        * run.length_scale**3
        * run.frequency
    )
    result = result.assign_coords(t=run.data.t / run.frequency)
    result.attrs = {
        "units": "W",
        "provenance": "reconstructed area-integrated outward energy faces",
    }
    return result


def compare_flux_diagnostics(run, fluxes):
    result = {}
    for name in fluxes:
        if "_flux_" in name or name not in run.data:
            continue
        error = fluxes[name] - run.data[name]
        if name == "ddt(Phi)" and fluxes.attrs["periodic_x"]:
            error -= error.isel(x=slice(run.mxg, -run.mxg)).mean(
                ("x", "z")
            )  # inversion gauge is free on each Y plane
        absolute = float(abs(run.interior(error)).max().compute())
        scale = float(abs(run.interior(run.data[name])).max().compute())
        result[name] = {
            "available": True,
            "max_absolute_error_normalised": absolute,
            "max_relative_error": absolute / max(scale, np.finfo(float).tiny),
            "max_integral_error_normalised": float(
                abs(run.integrate(error)).max().compute()
            ),
        }
    for name in (
        "DivEExB",
        "DivEReynolds",
        "DivEelectric",
        "DivEparallel",
        "DivEpol",
        "DivJpol",
        "phi_Div_Je",
        "electron_exb_work",
        "Eion_src",
        "ddt(Eion)",
    ):
        result.setdefault(name, {"available": False})
    return result
