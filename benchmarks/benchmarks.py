import dataclasses
import os
import pathlib
import shutil
import sys
import tempfile
from typing import ClassVar

from boututils.run_wrapper import launch_safe

REPO_BASE = pathlib.Path(__file__).parent.parent.resolve()
DEFAULT_TESTPATH = REPO_BASE / "tests" / "integrated"

sys.path.append(str(DEFAULT_TESTPATH))
from runtest_utils import AbstractDataDownloader, DataDownloader2DProduction


@dataclasses.dataclass
class BenchmarkParameters:
    testname: str
    nproc: int
    nout: int
    timestep: float
    testpath: pathlib.Path = DEFAULT_TESTPATH
    datadir: str = "data"
    data_downloader: type[AbstractDataDownloader] | None = None

    def __repr__(self) -> str:
        return self.testname


class Simulations:
    timeout = 300
    rounds = 3
    repeat = 1
    number = 1
    warmup_time = 0

    params: ClassVar = [
        BenchmarkParameters("1D-recycling-dthe", 1, 10, 2500),
        BenchmarkParameters(
            "2D-production", 10, 10, 10, data_downloader=DataDownloader2DProduction
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
        if param.data_downloader is not None:
            dd = param.data_downloader(param.testname, self.runpath)
            dd.download_data()
            dd.checkHash()
            dd.extractZip()

    def time_sim(self, param: BenchmarkParameters):
        launch_safe(
            f"hermes-3 -d {self.runpath / 'data'} nout={param.nout} timestep={param.timestep}",
            nproc=param.nproc,
        )

    def teardown(self, param: BenchmarkParameters):
        os.chdir(self.cwd)
        del self.rundir
