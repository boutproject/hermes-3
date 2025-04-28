from spack.package import *


class PyBoutdata(PythonPackage):
    """Python tools for working with BOUT++."""

    homepage = "https://github.com/boutproject/boutdata"
    pypi = "boutdata/boutdata-0.2.1.tar.gz"

    # maintainers("github_user1", "github_user2")

    # license("UNKNOWN", checked_by="github_user1")

    version(
        "0.2.1",
        sha256="043cddaeb38b128d2525f2005f48a5b7717ff5832a932183a4bef1d3eae389e0",
    )

    # Compatible Python versions
    depends_on("python@3.9:", type=("build", "run"))

    # Build dependencies
    depends_on("py-setuptools@61:", type="build")
    depends_on("py-setuptools-scm@6.2:+toml", type="build")

    # Runtime dependencies
    depends_on("py-boututils", type=("build", "run"))
    depends_on("py-matplotlib@3.2.1:", type=("build", "run"))
    depends_on("py-natsort@8.1.0:", type=("build", "run"))
    depends_on("py-netcdf4", type=("build", "run"))
    depends_on("py-numpy@1.22.0:", type=("build", "run"))
    depends_on("py-scipy@1.4.1:", type=("build", "run"))
    depends_on("py-sympy@1.5.1:", type=("build", "run"))

    def config_settings(self, spec, prefix):
        settings = {}
        return settings
