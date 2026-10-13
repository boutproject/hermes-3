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
    hermesconfig.config["src"] / "tests/integrated/1D-sheath-simple/__snapshots__"
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
        "Ee_sheath": norm(out["Ee_sheath"]),
        "Ei_sheath": norm(out["Ei_sheath"]),
        "Si_sheath": norm(out["Si_sheath"]),
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
        # Re-run with higher electron heat transmission.
        # Should give lower electron temperature and lower velocity.
        pytest.param("sheath_boundary_simple:gamma_e=35", id="Sheath higher gamma_e"),
        ################################################################################
        # Re-run with lower ion heat transmission.
        # Should give higher ion temperature and higher velocity.
        pytest.param("sheath_boundary_simple:gamma_i=0.35", id="Sheath lower gamma_i"),
        ################################################################################
        # Re-run with higher secondary electron emission.
        # Should give higher electron temperature and higher velocity.
        pytest.param(
            "sheath_boundary_simple:secondary_electron_coef=0.99", id="Sheath higher Ge"
        ),
        ################################################################################
        # Re-run with flipped direction. Should be mirror of base case
        pytest.param("lower_sheath=true", id="Mirrored BCs"),
        ################################################################################
        # Linear-free boundary
        pytest.param(
            "sheath_boundary_simple:density_boundary_mode=linear "
            "sheath_boundary_simple:pressure_boundary_mode=linear "
            "sheath_boundary_simple:temperature_boundary_mode=linear",
            id="Linear_free boundary mode",
        ),
    ],
)
def test_sheath_simple_case(extra_args, snaptolshot):
    """Run and check a whole test case"""

    snaptolshot.set_snaptolshot_dir(SNAPSHOT_DIR)

    _, out = launch_safe(f"{HERMES_BIN} -d {PATH} {extra_args}", nproc=1, pipe=True)

    Path("out.log").write_text(out)

    results = collect_all()

    assert snaptolshot(rtol=RTOL, atol=ATOL) == results


if __name__ == "__main__":
    import sys

    import pytest

    cli_args = sys.argv[1:] or ["-q"]

    pytest.main(args=[__file__] + cli_args)
