#!/usr/bin/env python3
"""Generate deterministic BOUT expressions for the continuum slab MMS."""

import argparse
import json
import tempfile
from pathlib import Path

import sympy as s
from boutdata.data import BoutOptionsFile
from boutdata.mms import exprToStr
from manufactured import (
    COORDINATES,
    EVOLVED,
    LENGTH,
    RECOVERED,
    VARIANTS,
    model,
    preflight,
)

HERE = Path(__file__).resolve().parent
REFERENCE_KEYS = dict(
    zip(
        RECOVERED,
        (
            "total_pressure",
            "pressure_i",
            "pressure_j",
            "ion_temperature",
            "potential",
            "pressure_potential",
            "current_potential",
        ),
    )
)


def input_path(variant):
    return HERE / f"data-{variant}" / "BOUT.inp"


def generate(destination, variant="baseline"):
    fields, _, sources, _ = model(variant)
    x, y, z, time = s.symbols("x y z t")
    replacement = dict(
        zip(
            COORDINATES,
            (LENGTH * x, LENGTH * y / (2 * s.pi), LENGTH * z / (2 * s.pi), time),
        )
    )
    keys = [(name, "solution") for name in EVOLVED]
    keys += [(name, "source") for name in EVOLVED]
    keys += [("mms", REFERENCE_KEYS[name]) for name in RECOVERED]
    expressions = [fields[name] for name in EVOLVED]
    expressions += [sources[name] for name in EVOLVED]
    expressions += [fields[name] for name in RECOVERED]
    expressions = [e.subs(replacement, simultaneous=True) for e in expressions]
    factors, reduced = s.cse(
        expressions, symbols=s.numbered_symbols("mms:c"), order="canonical"
    )
    # Context copies have a cost too. Inline cheap arithmetic and bind only
    # substantial repeated expressions (densities, gradients, pressure ratios).
    inline, substantial = {}, []
    for symbol, expression in factors:
        expression = expression.xreplace(inline)
        if s.count_ops(expression) < 12:
            inline[symbol] = expression
        else:
            substantial.append((symbol, expression))
    factors = substantial
    reduced = [expression.xreplace(inline) for expression in reduced]
    options = BoutOptionsFile(str(HERE / "base.inp"))
    options.getSection("mms")["variant"] = variant
    if variant == "skew":
        # The analytic chi and its normal derivative both vanish at the wall.
        # Dirichlet AC guards retain second-order mixed-metric wall transport;
        # plain X-Neumann guards leave a numerical g13*d_z(chi) wall current.
        options["evolve_vorticity_ion_energy"]["polarisation_ac_boundary"] = "dirichlet"
        for key, value in {
            "g13": "1/4",
            "g33": "17/16",
            "g_11": "17/16",
            "g_13": "-1/4",
        }.items():
            options["mesh"][key] = value
    if variant == "viscosity":
        # Units used by Hermes: nu_SI = nu_normalised * rho_s0^2 * Omega_ci.
        # With Tnorm=Bnorm=1, rho_s0^2*Omega_ci = Tnorm/Bnorm.
        options["evolve_vorticity_ion_energy"]["viscosity"] = "0.02"

    def exported(expression):
        result = exprToStr(expression)
        if (
            expression.has(s.Integral, s.Derivative, s.I)
            or any(symbol.name in ("X", "Y", "Z") for symbol in expression.free_symbols)
            or "**" in result
        ):
            raise ValueError(f"Unexportable MMS expression: {result}")
        return result

    # Use ordinary option references: BoutOptionsFile currently mistakes square
    # brackets in expression contexts for section headers when reading settings.
    for symbol, expression in factors:
        options.getSection("mms")[symbol.name.split(":")[1]] = exported(expression)

    for (section, key), expression in zip(keys, reduced):
        options.getSection(section)[key] = exported(expression)
        if section == "Vort" and key == "source":
            # Radial solver Vort guards carry fixed wall potentials. BOUT++
            # adds MMS forcing after the model RHS, including those guards.
            options[section][key] = f"where(x*(1-x),{exported(expression)},0)"
    for name in EVOLVED:
        if name == "Vort":
            continue  # Physical radial vorticity guards belong to radial_boundary
        boundary = (
            f"dirichlet_o2({name}:solution)" if name in ("Vort", "Eion") else "neumann"
        )
        for face in ("xin", "xout"):
            options[name]["bndry_" + face] = boundary
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    options.write(str(destination), overwrite=True)
    return {"expressions": len(expressions), "common_factors": len(factors)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--variant", choices=VARIANTS, default="baseline")
    parser.add_argument(
        "--check-files",
        action="store_true",
        help="Compare regenerated inputs with committed files",
    )
    parser.add_argument("--all", action="store_true", help="Generate every variant")
    parser.add_argument(
        "--check", action="store_true", help="Run analytic preflight checks"
    )
    args = parser.parse_args()
    if args.all and args.output:
        parser.error("--output cannot be combined with --all")
    reports = {}
    for variant in VARIANTS if args.all else (args.variant,):
        report = preflight(variant) if args.check else {}
        if args.check_files:
            with tempfile.TemporaryDirectory() as folder:
                regenerated = Path(folder) / "BOUT.inp"
                report.update(generate(regenerated, variant))
                if regenerated.read_text() != input_path(variant).read_text():
                    raise SystemExit(f"Stale MMS input: {input_path(variant)}")
                report["committed_input"] = "current"
        else:
            report.update(generate(args.output or input_path(variant), variant))
        reports[variant] = report
    print(json.dumps(reports, indent=2))


if __name__ == "__main__":
    main()
