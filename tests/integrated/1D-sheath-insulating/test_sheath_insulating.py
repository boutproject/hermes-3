#!/usr/bin/env python3
import os
import sys
from pathlib import Path

import pytest
from boutdata.collect import BoutArray
from boutdata.data import BoutOutputs
from boututils.run_wrapper import launch_safe

sys.path.append("../../../")

import hermesconfig

HERMES_BIN = os.environ.get("HERMES_BIN", hermesconfig.config["exe"])
PATH = "data"
SNAPSHOT_DIR = (
    hermesconfig.config["src"] / "tests/integrated/1D-sheath-insulating/__snapshots__"
)

# Loose tolerance due to low-res runs
RTOL = 1e-2
ATOL = 1e-2


def norm(var: BoutArray) -> BoutArray:
    """Normalise a variable by its conversion factor"""
    return var * var.attributes["conversion"]


def collect_all():
    """Collect all relevant variables from an output file"""
    out = BoutOutputs(PATH, tind=-1, caching=True)

    Pe = out["Pe"]
    Ne = out["Ne"]
    Tnorm = out["Tnorm"]
    Ti = out["Ti"] * Tnorm

    Te = (Pe / Ne) * Tnorm
    Te.attributes["long_name"] = "e temperature"
    Te.attributes["units"] = "eV"

    results = {
        "Te": Te,
        "Ti": Ti,
        "Ve": norm(out["Ve"]),
        "Vi": norm(out["Vi"]),
    }

    boundary_results = {"upstream": {}, "downstream": {}}

    for field, value in results.items():
        # Slices for upstream and downstream boundaries (interior points)
        boundary_results["upstream"][field] = value[(-1, 0, 0, 0)]
        boundary_results["downstream"][field] = value[(-1, 0, -1, 0)]

    return boundary_results


@pytest.mark.parametrize(
    "extra_args",
    [
        pytest.param("", id="Base case"),
        ################################################################################
        # Re-run with flipped direction. Should be mirror of base case
        pytest.param("lower_sheath=true", id="Mirrored BCs"),
    ],
)
@pytest.mark.snaptol_dir(SNAPSHOT_DIR)
def test_sheath_insulating(extra_args, snaptolshot):
    _, out = launch_safe(f"{HERMES_BIN} -d {PATH} {extra_args}", nproc=1, pipe=True)

    Path("out.log").write_text(out)

    results = collect_all()

    assert snaptolshot(rtol=RTOL, atol=ATOL) == results


if __name__ == "__main__":
    import sys

    import pytest

    cli_args = sys.argv[1:] or ["-q"]

    sys.exit(pytest.main(args=[__file__] + cli_args))
