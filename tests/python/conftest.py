import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))


import pytest


@pytest.fixture
def hermes_executable():
    root = Path(__file__).resolve().parents[2]
    executable = Path(
        os.environ.get("HERMES_EXECUTABLE", root / "build/hermes-3")
    ).resolve()
    if not executable.is_file():
        pytest.skip("Build Hermes-3 or set HERMES_EXECUTABLE to run native regressions")
    return executable


def pytest_addoption(parser):
    parser.addoption(
        "--run-slow",
        action="store_true",
        help="Run long native drift-plane regressions",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: long native integration")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--run-slow"):
        import pytest

        skip = pytest.mark.skip(reason="Use --run-slow for long native regressions")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip)
