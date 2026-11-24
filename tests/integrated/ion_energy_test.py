"""Native test launching; callers read and analyse the saved outputs separately."""

import argparse
import shlex
import shutil
import tempfile
from pathlib import Path

from boututils.run_wrapper import launch_safe


def parser(description):
    cli = argparse.ArgumentParser(description=description)
    root = Path(__file__).resolve().parents[2]
    cli.add_argument("--exe", type=Path, default=root / "build/hermes-3")
    cli.add_argument("--nproc", type=int, default=2)
    cli.add_argument(
        "--output-root",
        type=Path,
        default=Path(tempfile.mkdtemp(prefix="hermes-ion-energy-")),
    )
    return cli


def launch(args, label, template, overrides=(), *, nproc=None):
    directory = args.output_root / label
    directory.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(template, directory / "BOUT.inp")
    command = [
        str(args.exe.resolve()),
        "-d",
        str(directory.resolve()),
        "restart=false",
        *overrides,
    ]
    launch_safe(
        shlex.join(command),
        nproc=nproc or args.nproc,
        mthread=1,
        output=shlex.quote(str(directory / "run.log")),
    )
    return directory
