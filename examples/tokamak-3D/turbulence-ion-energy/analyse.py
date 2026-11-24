#!/usr/bin/env python3
"""Save power/particle tables and figures without launching a simulation."""

import argparse
from pathlib import Path

import xhermes

cli = argparse.ArgumentParser(description=__doc__)
cli.add_argument("directory", type=Path)
cli.add_argument("--output", type=Path)
cli.add_argument("--derivative", choices=("auto", "difference", "rhs"), default="auto")
args = cli.parse_args()
output = args.output or args.directory / "balance"
output.mkdir(parents=True, exist_ok=True)
with xhermes.open(args.directory, keep_xboundaries=True, keep_yboundaries=False) as ds:
    for quantity in (
        "energy",
        *["particles_" + name for name in ds.metadata["ion_species"]],
    ):
        result = ds.hermes.balance(quantity, derivative=args.derivative).compute()
        frame = xhermes.table(result)
        frame.to_csv(output / f"{quantity}.csv")
        result.to_dataframe().to_csv(output / f"{quantity}-history.csv")
        xhermes.plot_balance(result).savefig(output / f"{quantity}.pdf")
        print(f"{quantity}: {result.attrs['coverage']}, {result.attrs['derivative']}")
        if result.attrs["missing_terms"]:
            print("Missing: " + ", ".join(result.attrs["missing_terms"]))
        print(frame.to_string(float_format=lambda x: f"{x:.6g}"))
