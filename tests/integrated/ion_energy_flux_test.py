"""Shared implementation for independent stage-three source/flux test directories."""

import argparse
import json
import sys
import tempfile
from itertools import pairwise
from pathlib import Path

import numpy as np
import xarray as xr
from boutdata.data import BoutOptionsFile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from energy_analysis import (
    cadence_comparison,
    energy_balance,
    energy_inventory,
    integrated_power,
    open_run,
)
from energy_analysis.fluxes import (
    PrescribedSources,
    compare_flux_diagnostics,
    flux_power,
    net_face_flux,
    reconstruct_fluxes,
    reconstruct_sources,
    species_velocity,
)
from energy_analysis.report import inventory_report
from energy_analysis.transport import central_parallel
from ion_energy_balance import execute, histories, require, write_json

ENERGY_FAMILIES = ("DivEExB", "DivEReynolds", "DivEelectric", "DivEparallel", "DivEpol")


def validate_fixture(path, manifest):
    """Close source coverage from controlled input and an explicit SI manifest."""
    options = BoutOptionsFile(str(path))
    expected = manifest["case"]["components"]
    require(
        {c.strip() for c in options["hermes:components"].split(",")}
        == {"evolve_vorticity_ion_energy", *expected},
        "Unaccounted fixture component",
    )
    for name, types in expected.items():
        require(
            options[name]["type"] == types, f"Unaccounted fixture physics for {name}"
        )
        for component, prefix in (("evolve_density", "N"), ("evolve_pressure", "P")):
            if component in types.split(", "):
                value = options[prefix + name].evaluate_scalar("source")
                expected_rate = (
                    manifest["sources"]["thermal"][name] / 1.5 if prefix == "P" else 0
                )
                require(value == expected_rate, f"Source mismatch in {prefix + name}")
                for key in ("source_time_dependent", "source_only_in_core"):
                    require(
                        key not in options[prefix + name]
                        or not options[prefix + name].get_bool(key),
                        f"Unaccounted {key}",
                    )
        if "test_energy_sources" in types.split(", "):
            for kind, key in (
                ("thermal", "energy"),
                ("density", "density"),
                ("momentum", "momentum"),
            ):
                require(
                    options[name].evaluate_scalar("test_source_" + key)
                    == manifest["sources"][kind][name],
                    f"Prescribed {kind} source mismatch",
                )


def audit(run, manifest, diagnose):
    report, inventory = inventory_report(run)
    tol = manifest["tolerances"]
    sources = PrescribedSources(**manifest["sources"])
    flux = reconstruct_fluxes(run, sources=sources).compute()
    external = reconstruct_sources(run, sources)
    comparisons = compare_flux_diagnostics(run, flux)
    assertions = {
        "partition": report["partition_max_absolute_error_normalised"]
        < tol["partition_absolute"],
        "forward": report["forward_max_absolute_error_normalised"]
        < tol["forward_absolute"],
    }
    if diagnose:
        for name in (
            *ENERGY_FAMILIES,
            "DivJpol",
            "phi_Div_Je",
            "electron_exb_work",
            "Eion_src",
            "ddt(Eion)",
            "ddt(Pe)",
            "ddt(Phi)",
        ):
            check = comparisons.get(name, {"available": False})
            assertions["diagnostic_" + name] = (
                check["available"]
                and check["max_absolute_error_normalised"] < tol["diagnostic_absolute"]
            )
        for s in run.species:
            for component, prefix in (
                ("evolve_density", "N"),
                ("evolve_momentum", "NV"),
            ):
                if component in s.components:
                    name = "ddt(" + prefix + s.name + ")"
                    check = comparisons[name]
                    assertions["diagnostic_" + name] = (
                        check["max_absolute_error_normalised"]
                        < tol["diagnostic_absolute"]
                    )
        identity = (
            run.data["ddt(Eion)"]
            + sum(run.data[n] for n in ENERGY_FAMILIES)
            + run.data.phi_Div_Je
            - run.data.Eion_src
        )
        assertions["saved_ion_rhs_identity"] = (
            float(abs(run.interior(identity)).max()) < tol["diagnostic_absolute"]
        )
    gauss = {}
    for name in (*ENERGY_FAMILIES, "DivJpol"):
        difference = run.integrate(flux[name]) - net_face_flux(run, flux, name)
        error = float(abs(difference).max())
        gauss[name] = error
        assertions["gauss_" + name] = error < tol["gauss_absolute"]
    charge = sum(
        s.charge * flux["particle_pol_" + s.name] for s in run.species if s.name != "e"
    )
    assertions["species_charge_divergence"] = (
        float(abs(run.interior(charge - flux.DivJpol)).max())
        < tol["diagnostic_absolute"]
    )
    for dim in ("x", "z"):
        common = flux["DivJpol_flux_" + dim]
        split = sum(
            s.charge * flux["particle_pol_" + s.name + "_flux_" + dim]
            for s in run.species
            if s.name != "e"
        )
        assertions["species_charge_faces_" + dim] = (
            float(abs(run.interior(common - split)).max()) < tol["diagnostic_absolute"]
        )
    qion = integrated_power(run, flux.Eion_src - flux.phi_Div_Je)
    pion = sum(flux_power(run, flux, name) for name in ENERGY_FAMILIES).assign_attrs(
        units="W", provenance="reconstructed outward ion faces"
    )
    ion = energy_balance(
        inventory, reservoir="ion_saved", source_power=qion, outward_power=pion
    )
    # Full plasma: use external source rates; thermal/kinetic exchanges are
    # retained separately. The electron-inertia fixture has zero X-Z transport.
    outward = (
        pion
        + flux_power(run, flux, "electron_ExB")
        + flux_power(run, flux, "electron_parallel_thermal")
    )
    electron = next(s for s in run.species if s.name == "e")
    if "evolve_momentum" in electron.components:
        require(flux.attrs["periodic_y"], "Open electron-inertia budget is unavailable")
        require(
            float(abs(run.interior(run.data.phi - run.data.phi.isel(x=run.mxg))).max())
            < 1e-12,
            "Electron kinetic ExB flux is outside this test's coverage",
        )
        outward += flux_power(run, flux, "electron_parallel_kinetic")
    outward = outward.assign_attrs(
        units="W", provenance="reconstructed ion and electron faces"
    )
    plasma = energy_balance(
        inventory,
        source_power=integrated_power(run, external.plasma_src),
        outward_power=outward,
    )
    ve = species_velocity(run, electron)
    thermal_source = (
        external.thermal_e
        + flux.electron_exb_work
        + ve * central_parallel(run.data.Pe, run.data)
        + flux.electron_numerical_heating
    )
    electron_budget = energy_balance(
        inventory,
        reservoir="electron_thermal",
        source_power=integrated_power(
            run, thermal_source.assign_attrs(units_type="hermes")
        ),
        outward_power=(
            flux_power(run, flux, "electron_ExB")
            + flux_power(run, flux, "electron_parallel_thermal")
        ).assign_attrs(units="W"),
    )
    assertions["ion_cumulative"] = (
        float(abs(ion.relative_residual).max()) < tol["balance_relative"]
    )
    assertions["plasma_cumulative"] = (
        float(abs(plasma.relative_residual).max()) < tol["balance_relative"]
    )
    assertions["ion_intervals"] = (
        float(abs(ion.interval_relative_residual).max()) < tol["balance_relative"]
    )
    signals = {}
    for name, minimum in manifest["minimum_signals"].items():
        signals[name] = float(abs(run.interior(flux[name])).max())
        assertions["nonzero_" + name] = signals[name] > minimum
    if manifest["case"]["name"] == "parallel":
        # Preserve the uncorrected external-source/face balance. Quantify the
        # known BOUT heating defect separately; never classify it as heating input.
        defect = integrated_power(run, flux.electron_heating_defect)
        predicted = energy_balance(
            inventory,
            source_power=(
                integrated_power(run, external.plasma_src) + defect
            ).assign_attrs(
                units="W",
                provenance="external input plus tracked numerical heating defect",
            ),
            outward_power=outward,
        )
        del assertions["plasma_cumulative"]
        assertions["tracked_heating_defect_balance"] = (
            float(abs(predicted.relative_residual).max())
            < tol["heating_defect_relative"]
        )
        assertions["tracked_heating_defect_resolved"] = float(
            abs(plasma.relative_residual).max()
        ) > 20 * max(float(abs(predicted.relative_residual).max()), 1e-12)
        report["tracked_findings"] = {
            "parallel_numerical_heating_wave_coefficient": {
                "description": "Div_par_fvv_heating uses a +/- v_edge/2 instead of (a +/- v_edge)/2",
                "excess_power_W": defect.values.tolist(),
                "balance_with_predicted_defect": histories(predicted),
                "external_plasma_balance_conserves": False,
            }
        }
        thermal_work = integrated_power(
            run,
            (ve * central_parallel(run.data.Pe, run.data)).assign_attrs(
                units_type="hermes"
            ),
        )
        kinetic_work = integrated_power(run, flux.electron_kinetic_rhs)
        assertions["electron_kinetic_nonzero"] = (
            float(inventory.electron_parallel.max() / inventory.plasma_saved.max())
            > 1e-5
        )
        assertions["pressure_work_nonzero"] = float(abs(thermal_work).max()) > 1e-10
        report["electron_pressure_work_W"] = thermal_work.values.tolist()
        report["electron_kinetic_rhs_W"] = kinetic_work.values.tolist()
    if manifest["case"]["name"] == "heating":
        injected = float(abs(plasma.source_energy).max())
        assertions["deposition_resolved"] = injected > 100 * max(
            float(abs(plasma.residual).max()), 1e-16
        )
        expected = inventory.electron_thermal.isel(t=0) + plasma.source_energy
        assertions["analytic_thermal_deposition"] = (
            float(
                abs(inventory.electron_thermal - expected).max()
                / inventory.plasma_saved.max()
            )
            < tol["balance_relative"]
        )
    if manifest["case"]["name"] == "boundary":
        assertions["outward_power_resolved"] = float(
            plasma.outward_energy.isel(t=-1)
        ) > 100 * max(float(abs(plasma.residual).max()), 1e-16)
    if manifest["case"]["name"] == "sources":
        t = run.data.t / run.frequency
        expected = inventory.ion_parallel - inventory.ion_parallel.isel(t=0)
        deposited = (
            sum(sources.thermal[s.name] for s in run.species if s.name != "e")
            * float(run.integrate(xr.ones_like(run.data.Eion)).isel(t=0))
            * run.length_scale**3
        )
        expected += deposited * t.values
        assertions["analytic_source_energy"] = (
            float(
                abs(
                    inventory.ion_saved - inventory.ion_saved.isel(t=0) - expected
                ).max()
                / inventory.ion_saved.max()
            )
            < tol["balance_relative"]
        )
        for s in run.species:
            if "evolve_density" in s.components:
                expected_n = (
                    run.data["N" + s.name].isel(t=0)
                    + t * sources.density[s.name] / run.data.metadata["Nnorm"]
                )
                assertions["analytic_density_" + s.name] = (
                    float(abs(run.interior(expected_n - run.data["N" + s.name])).max())
                    < tol["history_absolute"]
                )
        for sp in run.species:
            if "evolve_momentum" in sp.components:
                saved = run.data["NV" + sp.name]
                expected_nv = (
                    saved.isel(t=0)
                    + t * sources.momentum[sp.name] / saved.attrs["conversion"]
                )
                assertions["analytic_momentum_" + sp.name] = (
                    float(abs(run.interior(saved - expected_nv)).max())
                    < tol["history_absolute"]
                )
    # Two disjoint X subvolumes partition the whole physical volume. Check
    # Gauss and transfer cancellation with independently selected cell sums.
    middle = run.data.sizes["x"] // 2
    regions = [slice(run.mxg, middle), slice(middle, run.data.sizes["x"] - run.mxg)]
    region_reports = []
    volume = run.data.J * run.data.dx * run.data.dy * run.data.dz
    for region in regions:
        select = {
            "x": region,
            "y": slice(run.myg, run.data.sizes["y"] - run.myg),
            "z": slice(0, run.data.sizes["z"]),
        }

        def integrate_region(value, select=select):
            return (value * volume).isel(select).sum(("x", "y", "z"), skipna=False)

        energy = (
            integrate_region(run.data.Eion) * run.energy_scale * run.length_scale**3
        )
        region_inventory = xr.Dataset({"ion_saved": energy}).assign_coords(
            t=inventory.t
        )
        region_inventory.ion_saved.attrs["units"] = "J"
        power = (
            (
                integrate_region(flux.Eion_src - flux.phi_Div_Je)
                * run.energy_scale
                * run.length_scale**3
                * run.frequency
            )
            .assign_coords(t=inventory.t)
            .assign_attrs(units="W")
        )
        boundary = sum(
            flux_power(run, flux, n, selection=select) for n in ENERGY_FAMILIES
        ).assign_attrs(units="W")
        budget = energy_balance(
            region_inventory,
            reservoir="ion_saved",
            source_power=power,
            outward_power=boundary,
        )
        assertions["region_" + str(region.start)] = (
            float(abs(budget.relative_residual).max())
            < tol["regional_balance_relative"]
        )
        region_reports.append(histories(budget))
    assertions["regional_face_cancellation"] = np.allclose(
        np.array(region_reports[0]["outward_energy"])
        + region_reports[1]["outward_energy"],
        ion.outward_energy,
        rtol=1e-10,
        atol=1e-16,
    )
    report.update(
        assertions=assertions,
        passed=all(assertions.values()),
        flux_coverage=flux.attrs,
        flux_diagnostics=comparisons,
        gauss_errors_normalised=gauss,
        signals_normalised=signals,
        ion_balance=histories(ion),
        plasma_balance=histories(plasma),
        electron_thermal_balance=histories(electron_budget),
        region_ion_balances=region_reports,
        source_power_W=integrated_power(run, external.plasma_src).values.tolist(),
        outward_power_W=outward.values.tolist(),
    )
    return (
        report,
        inventory,
        electron_budget if manifest["case"]["name"] == "parallel" else ion,
        flux,
    )


def main(inputs):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--profile", choices=("smoke", "ci"), default="ci")
    args = parser.parse_args()
    executable = args.exe or ROOT / "build/hermes-3"
    require(executable.is_file(), "Specify the built executable with --exe")
    manifest = json.loads((inputs / "criteria.json").read_text())
    fixture = inputs / "data/BOUT.inp"
    validate_fixture(fixture, manifest)
    output = Path(tempfile.mkdtemp(prefix=inputs.name + "-", dir=args.output_root))
    print(f"Results: {output}", flush=True)
    summary = {
        "manifest": manifest,
        "passed": False,
        "runs": [],
        "profile": args.profile,
    }
    write_json(output / "result.json", summary)
    cadence = manifest["cadences_normalised"]
    if args.profile == "smoke":
        cadence = cadence[-2:]
    runs, inventories, budgets, fluxes = {}, {}, {}, {}
    for dt, diagnose in [*((dt, True) for dt in cadence), (cadence[-1], False)]:
        label = f"dt{dt}-diagnose-{'on' if diagnose else 'off'}"
        folder = output / label
        print(f"Running {label}", flush=True)
        execute(
            executable.resolve(),
            folder,
            fixture,
            manifest["case"],
            dt,
            manifest["duration_normalised"],
            diagnose,
            args.timeout,
        )
        run = open_run(folder)
        run.data = run.data.load()
        report, inventory, budget, flux = audit(run, manifest, diagnose)
        write_json(folder / "energy-report.json", report)
        summary["runs"].append(str(folder))
        write_json(output / "result.json", summary)
        require(
            report["passed"],
            f"Failed {[n for n, v in report['assertions'].items() if not v]}: {folder / 'energy-report.json'}",
        )
        (
            runs[dt, diagnose],
            inventories[dt, diagnose],
            budgets[dt, diagnose],
            fluxes[dt, diagnose],
        ) = run, inventory, budget, flux
    fine = cadence[-1]
    xr.testing.assert_allclose(
        inventories[fine, True],
        inventories[fine, False],
        rtol=manifest["tolerances"]["history_relative"],
        atol=1e-13,
    )
    enabled, disabled = runs[fine, True], runs[fine, False]
    essential = {
        "Pe",
        "Eion",
        "Tion",
        "phi",
        "Vort",
        *["N" + s.name for s in enabled.species],
        *["NV" + s.name for s in enabled.species if "NV" + s.name in enabled.data],
        *["V" + s.name for s in enabled.species if "fixed_velocity" in s.components],
    }
    for name in essential:
        xr.testing.assert_allclose(
            enabled.data[name],
            disabled.data[name],
            rtol=manifest["tolerances"]["history_relative"],
            atol=manifest["tolerances"]["history_absolute"],
        )
    stripped = enabled.data.drop_vars(
        [
            name
            for name in enabled.data.data_vars
            if "t" in enabled.data[name].dims and name not in essential
        ]
    )
    stripped_flux = reconstruct_fluxes(
        open_run(stripped), sources=PrescribedSources(**manifest["sources"])
    )
    xr.testing.assert_allclose(stripped_flux, fluxes[fine, True])
    xr.testing.assert_allclose(
        energy_inventory(open_run(stripped)), inventories[fine, True]
    )
    refinement = [
        cadence_comparison(budgets[c, True], budgets[f, True])
        for c, f in pairwise(cadence)
    ]
    cadence_checks = {}
    reference = max(float(abs(inventories[fine, True].plasma_saved).max()), 1e-30)
    for i, comparison in enumerate(refinement):
        cadence_checks[f"sampling_error_{i}"] = (
            comparison["estimated_fine_sampling_error_J"] / reference
            < manifest["tolerances"]["maximum_sampling_relative"]
        )
    if manifest["case"]["name"] in ("fluxes", "sources") and len(refinement) > 1:
        cadence_checks["sampling_refinement"] = (
            refinement[-1]["max_residual_difference_J"]
            < manifest["tolerances"]["maximum_sampling_refinement_ratio"]
            * refinement[-2]["max_residual_difference_J"]
        )
    summary["cadence_assertions"] = cadence_checks
    write_json(output / "result.json", summary)
    require(
        all(cadence_checks.values()), f"Cadence audit failed: {output / 'result.json'}"
    )
    summary.update(
        cadence_refinement=refinement,
        diagnostics_independent=True,
        optional_fields_stripped=True,
        passed=True,
    )
    write_json(output / "result.json", summary)
    print(f"{inputs.name} passed (flux/source audits and diagnostics pairs)")
    return 0


def run_test(inputs):
    try:
        return main(inputs)
    except (ValueError, AssertionError) as error:
        print(f"Energy flux integration test failed: {error}", file=sys.stderr)
        return 1
