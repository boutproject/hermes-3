import dataclasses
import os
import pathlib
import shutil
import sys
import tempfile
from collections.abc import Callable
from typing import ClassVar

from boututils.run_wrapper import launch_safe

REPO_BASE = pathlib.Path(__file__).parent.parent.resolve()


@dataclasses.dataclass
class BenchmarkParameters:
    testname: str
    nproc: int
    nout: int
    timestep: float
    testpath: pathlib.Path = REPO_BASE / "tests" / "integrated"
    datadir: str = "data"
    get_data: Callable[["BenchmarkParameters", pathlib.Path], None] = (
        lambda self, path: None
    )

    def __repr__(self) -> str:
        return self.testname


def get_2d_production_data(self: BenchmarkParameters, runpath: pathlib.Path):
    sys.path.append(str(self.testpath / self.testname))
    import runtest_utils

    sys.path.pop()
    dd = runtest_utils.DataDownloader(runpath)
    dd.download_data()
    dd.checkHash()
    dd.extractZip()


class TimingBenchmark:
    timeout = 180
    rounds = 3
    repeat = 1
    number = 1
    warmup_time = 0

    params: ClassVar = [
        BenchmarkParameters("1D-recycling-dthe", 1, 10, 2500),
        BenchmarkParameters(
            "2D-production", 10, 10, 10, get_data=get_2d_production_data
        ),
        BenchmarkParameters(
            "tokamak-2D-driftplane",
            4,
            15,
            250,
            testpath=REPO_BASE / "examples",
            datadir="2D-drift-plane-turbulence-te-ti",
        ),
    ]
    param_names: ClassVar = ["testcase"]

    def setup(self, param: BenchmarkParameters):
        self.rundir = tempfile.TemporaryDirectory()
        self.runpath = pathlib.Path(self.rundir.name)
        shutil.copytree(
            param.testpath / param.testname / param.datadir, self.runpath / "data"
        )
        self.cwd = os.getcwd()
        os.chdir(self.runpath)
        param.get_data(param, self.runpath)

    def time_sim(self, param: BenchmarkParameters):
        launch_safe(
            f"hermes-3 -d {self.runpath / 'data'} nout={param.nout} timestep={param.timestep}",
            nproc=param.nproc,
        )

    def teardown(self, param: BenchmarkParameters):
        os.chdir(self.cwd)
        del self.rundir
