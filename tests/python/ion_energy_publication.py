"""Run the fixed-wall D=0.1 publication input with its matched startup state."""

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "examples/tokamak-2D-driftplane/2D-energy/publication_scan.py"
spec = importlib.util.spec_from_file_location("ion_energy_publication_scan", SOURCE)
scan = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan)


def run_case(executable, folder, n, *, endpoint, output_step, nproc=2):
    """Use original tolerances and fixed-width pressure wall buffers."""
    outputs = round(endpoint / output_step)
    first_snapshot = endpoint / 4 if outputs % 4 == 0 else endpoint / 5
    config = scan.configuration(
        scan.parser().parse_args(
            [
                "--levels",
                "3",
                "--base-nx",
                "16",
                "--base-nz",
                "16",
                "--custom-D",
                "0.1",
                "--phi-boundary",
                "fixed",
                "--nproc",
                str(nproc),
                "--end-time",
                str(endpoint),
                "--outputs",
                str(outputs),
                "--convergence-time",
                str(endpoint),
                "--snapshot-times",
                "0",
                str(first_snapshot),
                str(endpoint * 0.5),
                str(endpoint),
            ]
        )
    )
    folder.mkdir(parents=True, exist_ok=True)
    row = {"formulation": "total_energy", "boundary": "fixed", "nx": n, "nz": n}
    (folder / "BOUT.inp").write_text(scan.input_text(config, row))
    (folder / "configuration.json").write_text(json.dumps(config, indent=2) + "\n")
    scan.build_startup(executable, folder, config, timeout=1200)
    scan.execute(executable, folder, config, timeout=7200)
    return folder
