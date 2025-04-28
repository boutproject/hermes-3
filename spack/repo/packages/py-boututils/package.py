from spack.package import *


class PyBoututils(PythonPackage):
    """pip-package of what was previously found in BOUT-dev/tools/pylib/boututils"""

    homepage = "https://www.example.com"
    pypi = "boututils/boututils-0.2.1.tar.gz"

    # maintainers("github_user1", "github_user2")

    # license("UNKNOWN", checked_by="github_user1")

    version(
        "0.2.1",
        sha256="b8e12cace0638645de09647b60f1f1e0079ded359855a9d220aede94736d2fe5",
    )

    # Compatible Python versions
    depends_on("python@3.8:", type=("build", "run"))

    # Build dependencies
    depends_on("py-setuptools@61:", type="build")
    depends_on("py-setuptools-scm@6.2:+toml", type="build")

    # Runtime dependencies
    depends_on("py-matplotlib@3.2.1:", type=("build", "link", "run"))
    depends_on("py-netcdf4", type=("build", "link", "run"))
    depends_on("py-numpy@1.22.0:", type=("build", "link", "run"))
    depends_on("py-scipy@1.4.1:", type=("build", "link", "run"))

    def config_settings(self, spec, prefix):
        settings = {}
        return settings
