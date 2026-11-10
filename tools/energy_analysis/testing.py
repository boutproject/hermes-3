"""Shared run provenance, launch, comparison and reporting for energy tests."""

import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np

from .dataset import AuditError, open_run
from .report import revision

# BOUT++'s output monitor prints these rows only after saving an output.
OUTPUT_TIME = re.compile(r"^\s*([+-]?\d+\.\d+[eE][+-]?\d+)\s+\d+\s+", re.MULTILINE)
FLOOR_MESSAGE = re.compile(r"\bfloor(?:ed|ing)?\b", re.IGNORECASE)


def run_progress(contents):
    """Read saved output times and runtime floor messages from a Hermes log.

    Configuration messages such as 'Option ... density_floor' are not floor
    events. Hermes currently clips density and recovered pressure silently;
    this also recognises explicit 'Ion energy below floor' diagnostics.
    """
    lines = contents.replace("\r", "\n")
    times = [float(match[1]) for match in OUTPUT_TIME.finditer(lines)]
    floor = next(
        (
            line.strip()
            for line in lines.splitlines()
            if not line.lstrip().startswith("Option ") and FLOOR_MESSAGE.search(line)
        ),
        None,
    )
    return max(times) if times else None, floor


def stop_run(process):
    """Reap the launcher and its MPI ranks when a monitored run fails."""
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def parser(description):
    result = argparse.ArgumentParser(description=description)
    root = Path(__file__).resolve().parents[2]
    result.add_argument("--exe", type=Path, default=root / "build/hermes-3")
    result.add_argument("--output-root", type=Path)
    result.add_argument("--timeout", type=float, default=180)
    result.add_argument("--nproc", type=int, default=1)
    result.add_argument("--mpiexec", default=shutil.which("mpiexec"))
    result.add_argument("--mpi-args", default="")
    return result


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def resolve_limit(criterion, run=None, *, minimum=False, truncation=False):
    """Resolve a numeric floor or a solver-relative limit from run options."""
    if minimum or truncation:
        limit = float(criterion)
    else:
        specification = (
            criterion if isinstance(criterion, dict) else {"floor": criterion}
        )
        if set(specification) - {"floor", "rtol_factor"}:
            raise AuditError(f"Unknown limit fields: {specification}")
        floor = float(specification["floor"])
        factor = float(specification.get("rtol_factor", 0))
        if not np.isfinite(floor) or not np.isfinite(factor) or factor < 0:
            raise AuditError(f"Invalid limit: {specification}")
        if factor:
            if run is None:
                raise AuditError("An rtol-scaled limit requires a run")
            try:
                rtol = float(run.data.options["solver"].evaluate_scalar("rtol"))
            except (KeyError, ValueError, TypeError, NameError, SyntaxError) as error:
                raise AuditError(f"Cannot resolve solver:rtol: {error}") from error
            if not np.isfinite(rtol) or rtol <= 0:
                raise AuditError(f"Invalid solver:rtol: {rtol}")
            limit = max(floor, factor * rtol)
        else:
            limit = floor
        if limit < 1e-10:
            raise AuditError(f"Effective limit {limit} is below 1e-10")
    if not np.isfinite(limit) or limit < 0:
        raise AuditError(f"Invalid threshold: {limit}")
    return limit


class Study:
    def __init__(self, here, args):
        self.here, self.args = Path(here), args
        if args.output_root:
            args.output_root.mkdir(parents=True, exist_ok=True)
        self.output = Path(
            tempfile.mkdtemp(prefix=self.here.name + "-", dir=args.output_root)
        )
        self.criteria = json.loads((self.here / "criteria.json").read_text())
        self.result = {
            "passed": False,
            "checks": [],
            "cases": {},
            "analysis": revision(Path(__file__).resolve().parents[2]),
        }
        self.start = time.monotonic()
        self.current_run = None
        print(f"Results: {self.output}", flush=True)

    def __enter__(self):
        return self

    def __exit__(self, kind, error, traceback):
        self.result.update(
            passed=kind is None
            and all(c["passed"] for c in self.result["checks"])
            and all(
                c.get("outcome") == "completed" for c in self.result["cases"].values()
            ),
            elapsed_seconds=time.monotonic() - self.start,
        )
        if error:
            self.result["error"] = str(error)
        write_json(self.output / "result.json", self.result)
        if kind is None and not self.result["passed"]:
            raise AuditError("Criteria failed; see " + str(self.output / "result.json"))

    def check(
        self,
        name,
        value,
        criterion=None,
        *,
        minimum=False,
        limit=None,
        run=None,
        truncation=False,
    ):
        if limit is None:
            limit = self.criteria[criterion or name]
        run = run if run is not None else self.current_run
        limit = resolve_limit(limit, run, minimum=minimum, truncation=truncation)
        value = float(value)
        passed = bool(
            np.isfinite(value) and (value >= limit if minimum else value <= limit)
        )
        self.result["checks"].append(
            {
                "name": name,
                "value": value,
                "limit": limit,
                "comparison": ">=" if minimum else "<=",
                "kind": "minimum"
                if minimum
                else "truncation"
                if truncation
                else "limit",
                "passed": passed,
            }
        )
        print(
            f"{'PASS' if passed else 'FAIL'} {name}: {value:.6g} {'>=' if minimum else '<='} {limit:.6g}",
            flush=True,
        )

    def run(
        self, name, overrides=(), *, template=None, restart=None, stall_seconds=None
    ):
        folder = self.output / name
        folder.mkdir()
        template = Path(template or self.here / "data/BOUT.inp")
        shutil.copyfile(template, folder / "BOUT.inp")
        if restart:
            for file in Path(restart).glob("BOUT.restart.*.nc"):
                shutil.copyfile(file, folder / file.name)
        self.launch(folder, overrides, stall_seconds=stall_seconds)
        self.current_run = open_run(folder)
        self.result["cases"][name]["solver_rtol"] = float(
            self.current_run.data.options["solver"].evaluate_scalar("rtol")
        )
        return self.current_run

    def launch(self, folder, overrides=(), *, stall_seconds=None):
        args = self.args
        stall_seconds = float(
            self.criteria.get("stall_seconds", 120)
            if stall_seconds is None
            else stall_seconds
        )
        if not np.isfinite(stall_seconds) or stall_seconds <= 0:
            raise AuditError(f"Invalid stall_seconds: {stall_seconds}")
        # BOUT++ rejects duplicate command-line keys. Later overrides replace
        # defaults, particularly restart=false in a restart round trip.
        overrides = list(
            {str(value).split("=", 1)[0]: str(value) for value in overrides}.values()
        )
        command = [
            str(args.exe.resolve()),
            "-d",
            str(folder.resolve()),
            *map(str, overrides),
        ]
        if args.nproc > 1:
            if not args.mpiexec:
                raise AuditError("--nproc requires an MPI launcher")
            command = [
                str(args.mpiexec),
                *shlex.split(args.mpi_args),
                "-n",
                str(args.nproc),
                *command,
            ]
        status = {
            "command": command,
            "input_sha256": hashlib.sha256(
                (folder / "BOUT.inp").read_bytes()
            ).hexdigest(),
            "executable_sha256": hashlib.sha256(args.exe.read_bytes()).hexdigest(),
            "stall_seconds": stall_seconds,
            "last_output_time": None,
            "maximum_output_interval_seconds": 0.0,
            "output_progress": [],
            "outcome": "failed",
        }
        start = time.monotonic()
        last_advance = start

        def record_progress(latest, now):
            nonlocal last_advance
            if latest is not None and (
                status["last_output_time"] is None
                or latest > status["last_output_time"]
            ):
                if status["last_output_time"] is not None:
                    status["maximum_output_interval_seconds"] = max(
                        status["maximum_output_interval_seconds"], now - last_advance
                    )
                status["last_output_time"] = latest
                status["output_progress"].append(
                    {"time": latest, "elapsed_seconds": now - start}
                )
                last_advance = now

        self.result["cases"][folder.name] = status
        log_path = folder / "run.log"
        with log_path.open("w") as log:
            try:
                process = subprocess.Popen(
                    command,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
                try:
                    while True:
                        latest, floor_message = run_progress(
                            log_path.read_text(errors="replace")
                        )
                        now = time.monotonic()
                        record_progress(latest, now)
                        if floor_message:
                            status.update(
                                outcome="floored", floor_message=floor_message
                            )
                            stop_run(process)
                            break
                        if process.poll() is not None:
                            # Read once more after exit, including the final flushed rows.
                            latest, floor_message = run_progress(
                                log_path.read_text(errors="replace")
                            )
                            record_progress(latest, time.monotonic())
                            if floor_message:
                                status.update(
                                    outcome="floored", floor_message=floor_message
                                )
                            elif (
                                process.returncode == 0
                                and len(list(folder.glob("BOUT.dmp.*.nc")))
                                == args.nproc
                            ):
                                status["outcome"] = "completed"
                            break
                        if now - last_advance >= stall_seconds:
                            status["outcome"] = "stalled"
                            stop_run(process)
                            break
                        if now - start >= args.timeout:
                            status["timed_out"] = True
                            stop_run(process)
                            break
                        time.sleep(min(1.0, stall_seconds / 10, args.timeout / 10))
                finally:
                    if process.poll() is None:
                        stop_run(process)
                status["returncode"] = process.returncode
            finally:
                status["elapsed_seconds"] = time.monotonic() - start
                write_json(folder / "execution.json", status)
                write_json(self.output / "result.json", self.result)
        if status["outcome"] != "completed":
            self.check("execution/" + folder.name, 1, limit=0, truncation=True)
            reason = "timed out" if status.get("timed_out") else status["outcome"]
            raise AuditError(
                f"Run {reason} (last output {status['last_output_time']}): {log_path}"
            )


def state_difference(first, second, *, first_index=-1, second_index=-1):
    names = ["Eion", "Vort", "phi", "Tion", "Pe"]
    names += ["N" + s.name for s in first.species]
    names += ["NV" + s.name for s in first.species if "NV" + s.name in first.data]
    return max(
        float(
            abs(
                first.data[n].isel(t=first_index, drop=True)
                - second.data[n].isel(t=second_index, drop=True)
            ).max()
        )
        for n in names
    )
