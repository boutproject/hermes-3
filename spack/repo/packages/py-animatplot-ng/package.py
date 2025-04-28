from spack.package import *


class PyAnimatplotNg(PythonPackage):
    """A python package for making interactive as well as animated plots with matplotlib. Based on animatplot by r-makaro."""

    homepage = "https://github.com/boutproject/animatplot-ng/"
    pypi = "animatplot-ng/animatplot-ng-0.4.4.tar.gz"

    # maintainers("github_user1", "github_user2")

    # license("UNKNOWN", checked_by="github_user1")

    version(
        "0.4.4",
        sha256="89f51ca4d63714a918f95ef14d576f420ae6f2aad08968e634379634ca375324",
    )

    # Compatible Python versions
    depends_on("python@3.5:", type=("build", "run"))

    # Build dependencies
    depends_on("py-setuptools@42:", type="build")
    depends_on("py-setuptools-scm@7:+toml", type="build")
    depends_on("py-wheel@0.29:", type="build")

    # Runtime dependencies
    depends_on("py-matplotlib@2.2:", type=("build", "run"))

    def config_settings(self, spec, prefix):
        settings = {}
        return settings
