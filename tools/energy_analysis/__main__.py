"""Print inventories and balances, including diagnostic-independent coverage."""

import argparse
import json
from pathlib import Path

from .dataset import AuditError, open_run
from .report import summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(argv)
    try:
        run = open_run(args.path)
        try:
            report = summary(run)
        finally:
            run.data.close()
    except (AuditError, OSError) as error:
        parser.exit(1, f"Energy audit failed: {error}\n")
    output = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.json:
        args.json.write_text(output)
    print(output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
