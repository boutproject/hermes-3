#!/usr/bin/env python3
"""Standalone runner copied into each prepared case as run-remote.py."""

import argparse
import hashlib
import json
import os
import platform
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def runtime_artifacts(exe):
    artifacts = {exe.resolve()}
    for directory in (
        exe.parent,
        exe.parent / "lib",
        exe.parent / "external/BOUT-dev/lib",
    ):
        for prefix in ("libhermes-3-lib", "libbout++"):
            for path in directory.glob(prefix + "*"):
                if path.is_file() and (
                    path.suffix in (".dylib", ".dll") or ".so" in path.name
                ):
                    artifacts.add(path.resolve())
    return {str(path): digest(path) for path in sorted(artifacts)}


def revision(path):
    def git(*args):
        result = subprocess.run(
            ["git", "-C", str(path), *args], capture_output=True, text=True, check=False
        )
        return result.stdout.strip() if result.returncode == 0 else None

    try:
        commit = git("rev-parse", "HEAD")
        dirty = git("status", "--porcelain", "--untracked-files=no")
        return {"revision": commit, "dirty": bool(dirty) if dirty is not None else None}
    except OSError:
        return {"revision": None, "dirty": None}


def main(argv=None, *, directory=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--mpiexec", default="mpiexec")
    p.add_argument(
        "--source-root",
        type=Path,
        help="Remote Hermes checkout for revision provenance",
    )
    p.add_argument("--timeout", type=float, help="Optional wall-clock limit in seconds")
    args = p.parse_args(argv)
    directory = Path(directory or Path(__file__).resolve().parent).resolve()
    case_file = directory / "remote-case.json"
    case = json.loads(case_file.read_text())
    if case["schema_version"] != 1:
        raise ValueError("Unsupported remote case schema")
    for name, expected in case["input_files"].items():
        if digest(directory / name) != expected:
            raise ValueError(f"Prepared input changed: {name}")
    if (
        list(directory.glob("BOUT.dmp.*.nc"))
        or (directory / "remote-execution.json").exists()
    ):
        raise ValueError("Outputs already exist; use a fresh copy of the prepared case")
    exe = args.exe.resolve()
    if not exe.is_file():
        raise ValueError(f"Executable not found: {exe}")
    if args.timeout is not None and args.timeout <= 0:
        raise ValueError("Timeout must be positive")
    command = [str(exe), "-d", str(directory)]
    if case["configuration"]["nproc"] > 1:
        command = [args.mpiexec, "-n", str(case["configuration"]["nproc"]), *command]
    artifacts = runtime_artifacts(exe)
    source = args.source_root.resolve() if args.source_root else None
    if source is None:
        try:
            result = subprocess.run(
                ["git", "-C", str(exe.parent), "rev-parse", "--show-toplevel"],
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode == 0:
                source = Path(result.stdout.strip())
        except OSError:
            pass
    status = {
        "schema_version": 1,
        "case": case["case"]["directory"],
        "prepared_case_sha256": digest(case_file),
        "input_files": case["input_files"],
        "command": command,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python": sys.version,
        "runtime_artifacts": artifacts,
        "revisions": {
            "hermes": revision(source) if source else {"revision": None, "dirty": None},
            "bout": revision(source / "external/BOUT-dev")
            if source
            else {"revision": None, "dirty": None},
        },
    }
    for name, actual in status["revisions"].items():
        expected = case.get("expected_revisions", {}).get(name, {}).get("revision")
        if expected and actual["revision"] and expected != actual["revision"]:
            raise ValueError(f"Remote {name} revision differs from the prepared study")
    start = time.monotonic()
    interrupted = False
    with (directory / "run.log").open("w") as stream:
        process = subprocess.Popen(
            command, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True
        )
        try:
            code = process.wait(timeout=args.timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
            interrupted = isinstance(error, KeyboardInterrupt)
            status["termination_reason"] = "interrupted" if interrupted else "timeout"
            code = 130 if interrupted else 124
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
    status.update(
        exit_code=code,
        elapsed_s=time.monotonic() - start,
        finished_utc=datetime.now(timezone.utc).isoformat(),
        runtime_artifacts_after=runtime_artifacts(exe),
    )
    path = directory / "remote-execution.json"
    temp = path.with_suffix(".json.tmp")
    temp.write_text(json.dumps(status, indent=2) + "\n")
    temp.replace(path)
    print(
        f"Exit code {code}; copy dumps, logs and remote-execution.json back for validation."
    )
    return 130 if interrupted else code


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError) as error:
        print(f"Remote run stopped: {error}", file=sys.stderr)
        sys.exit(1)
