"""Command line inventory audit; optional sources/fluxes are never assumed zero."""

import argparse
import json
from pathlib import Path

from .dataset import AuditError, open_run
from .fluxes import PrescribedSources
from .report import flux_report, inventory_report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", nargs="?", default="data-Eion")
    parser.add_argument("--json", type=Path, help="write a compact JSON report")
    parser.add_argument(
        "--fluxes", action="store_true", help="audit supported slab fluxes and sources"
    )
    parser.add_argument(
        "--sources",
        type=Path,
        help="JSON PrescribedSources manifest of SI rates and provenance",
    )
    args = parser.parse_args(argv)
    try:
        run = open_run(args.path)
        if args.fluxes or args.sources:
            sources = (
                PrescribedSources(**json.loads(args.sources.read_text()))
                if args.sources
                else None
            )
            report = flux_report(run, sources=sources)
        else:
            report, _ = inventory_report(run)
    except AuditError as error:
        parser.exit(1, f"Energy audit failed: {error}\n")
    output = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.json:
        args.json.write_text(output)
        print(f"Inventory report: {args.json}")
        print(report["balance"]["coverage"])
    else:
        print(output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
