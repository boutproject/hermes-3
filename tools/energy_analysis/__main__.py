"""Command line inventory audit; optional sources/fluxes are never assumed zero."""

import argparse
import json
from pathlib import Path

from .dataset import AuditError, open_run
from .report import inventory_report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", nargs="?", default="data-Eion")
    parser.add_argument("--json", type=Path, help="write a compact JSON report")
    args = parser.parse_args(argv)
    try:
        report, _ = inventory_report(open_run(args.path))
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
