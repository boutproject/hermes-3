"""Compact machine-readable inventory reports."""

import subprocess
import sys
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path

import numpy as np

from .balance import energy_balance, regular_timestep
from .inventory import compare_diagnostics, energy_densities, energy_inventory


def revision(path):
    def git(*args):
        result = subprocess.run(
            ["git", "-C", str(path), *args], capture_output=True, text=True, check=False
        )
        return result.stdout.strip() if result.returncode == 0 else None

    return {
        "revision": git("rev-parse", "HEAD"),
        "dirty": bool(git("status", "--porcelain")),
    }


def inventory_report(run):
    densities = energy_densities(run)
    inventory = energy_inventory(run, densities=densities).compute()
    dt = regular_timestep(inventory.t)
    partition = run.interior(densities.ion_reconstructed - run.data.Eion)
    forward = run.interior(densities.vorticity_reconstructed - run.data.Vort)
    import xhermes

    checkout = Path(xhermes.__file__).resolve().parents[1]
    report = {
        "schema_version": 1,
        "source": run.source,
        "python": sys.executable,
        "xhermes": {
            "version": version("xhermes"),
            "path": str(checkout),
            **revision(checkout),
        },
        "hermes_revision": run.data.metadata.get("HERMES_REVISION"),
        "analysis_checkout": revision(Path(__file__).resolve().parents[2]),
        "run_id": run.data.metadata.get("run_id"),
        "normalisations": {
            "energy_density_J_m3": run.energy_scale,
            "length_m": run.length_scale,
            "frequency_s_inv": run.frequency,
        },
        "domain": {"mxg": run.mxg, "myg": run.myg, "dimensions": dict(run.data.sizes)},
        "species": [asdict(s) for s in run.species],
        "time_s": inventory.t.values.tolist(),
        "time_spacing_s": dt,
        "inventory_J": {name: inventory[name].values.tolist() for name in inventory},
        "partition_max_absolute_error_normalised": float(
            abs(partition).max().compute()
        ),
        "partition_error_by_time_normalised": abs(partition)
        .max(["x", "y", "z"])
        .compute()
        .values.tolist(),
        "forward_max_absolute_error_normalised": float(abs(forward).max().compute()),
        "forward_error_by_time_normalised": abs(forward)
        .max(["x", "y", "z"])
        .compute()
        .values.tolist(),
        "diagnostics": compare_diagnostics(run, densities=densities),
        "balance": dict(energy_balance(inventory).attrs),
    }
    if not all(np.isfinite(inventory[name]).all() for name in inventory):
        from .dataset import AuditError

        raise AuditError("Nonfinite reconstructed inventory; check saved guard fields")
    return report, inventory
