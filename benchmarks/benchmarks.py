import os
import pathlib
import sys
import tempfile
import shutil

from boututils.run_wrapper import launch_safe

TESTPATH = pathlib.Path(__file__).parent.parent.resolve() / "tests" / "integrated"


class Time1DRun:
    """
    Simple performance benchmark for 1D case.
    """

    timeout = 180
    rounds = 1
    repeat = 1
    number = 1
    warmup_time = 0

    def setup(self):
        self.rundir = tempfile.TemporaryDirectory()
        self.runpath = pathlib.Path(self.rundir.name)
        # FIXME: is this the right test? My notes say 1D-thene
        shutil.copytree(TESTPATH / "1D-recycling-dthe" / "data", self.runpath / "data")
        self.cwd = os.getcwd()
        os.chdir(self.runpath)

    def time_1D_sim(self):
        launch_safe(
            f"hermes-3 -d {self.runpath / 'data'} nout=10 timestep=2500", nproc=1
        )

    def teardown(self):
        os.chdir(self.cwd)
        del self.rundir


class Time2DRun:
    """
    Simple performance benchmark based on running 2D-production
    """

    timeout = 180
    rounds = 1
    repeat = 1
    number = 1
    warmup_time = 0

    def setup(self):
        self.rundir = tempfile.TemporaryDirectory()
        self.runpath = pathlib.Path(self.rundir.name)
        testdir = TESTPATH / "2D-production"
        shutil.copytree(testdir / "data", self.runpath / "data")
        self.cwd = os.getcwd()
        os.chdir(self.runpath)
        # Load the runtest file to access methods needed to download data
        sys.path.append(str(testdir))
        # print(testdir / "runtest")
        # spec = importlib.util.spec_from_file_location("runtest", testdir / "runtest")
        # assert spec is not None
        # assert spec.loader is not None
        # runtest = importlib.util.module_from_spec(spec)
        # sys.modules["runtest"] = runtest
        # spec.loader.exec_module(runtest)
        # Download test data
        import runtest_utils

        dd = runtest_utils.DataDownloader(self.runpath)
        dd.download_data()
        dd.checkHash()
        dd.extractZip()

    def time_2D_sim(self):
        launch_safe(
            f"hermes-3 -d {self.runpath / 'data'} nout=10 timestep=10", nproc=10
        )

    def teardown(self):
        os.chdir(self.cwd)
        del self.rundir
