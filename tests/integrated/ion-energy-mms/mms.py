#!/usr/bin/env python3
"""Generate deterministic BOUT expressions for the continuum slab MMS."""

import argparse
import json
from pathlib import Path

import sympy as s
from boutdata.data import BoutOptionsFile
from boutdata.mms import exprToStr
from manufactured import COORDINATES, EVOLVED, LENGTH, RECOVERED, model, preflight

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


def generate(destination):
    fields, _, sources, _ = model()
    names = EVOLVED
    x, y, z, time = s.symbols("x y z t")
    replacement = dict(
        zip(
            COORDINATES,
            (LENGTH * x, LENGTH * y / (2 * s.pi), LENGTH * z / (2 * s.pi), time),
        )
    )
    keys = [(name, "solution") for name in names]
    keys += [(name, "source") for name in names]
    keys += [("mms", REFERENCE_KEYS[name]) for name in RECOVERED]
    expressions = [fields[name] for name in names]
    expressions += [sources[name] for name in names]
    expressions += [fields[name] for name in RECOVERED]
    keys += [
        ("mms", key)
        for key in (
            "electron_temperature",
            "electron_velocity",
            "parallel_current",
            "parallel_current_divergence",
        )
    ]
    expressions += [fields[key] for key in ("Te", "Ve", "Jpar", "DivJpar")]
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
    for name in names:
        if name == "Vort":
            continue  # Physical radial vorticity guards belong to radial_boundary
        boundary = (
            f"dirichlet_o2({name}:solution)" if name in ("Eion", "NVe") else "neumann"
        )
        for face in ("xin", "xout"):
            options[name]["bndry_" + face] = boundary
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    options.write(str(destination), overwrite=True)
    return {"expressions": len(expressions), "common_factors": len(factors)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "data/BOUT.inp")
    parser.add_argument(
        "--check", action="store_true", help="Check continuum identities"
    )
    args = parser.parse_args()
    report = preflight() if args.check else {}
    report.update(generate(args.output))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
