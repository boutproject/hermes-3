from spack.package import *


class PyGelidum(PythonPackage):
    """Inspired by the method freeze found in other languages like Javascript, this package tries to make immutable objects to make it easier avoiding accidental modifications in your code."""

    homepage = "https://github.com/diegojromerolopez/gelidum"
    pypi = "gelidum/gelidum-0.8.2.tar.gz"

    # maintainers("github_user1", "github_user2")

    # license("UNKNOWN", checked_by="github_user1")

    version(
        "0.8.2",
        sha256="3761191eeb11a406620bcbc853730bfa82b2b947cb55b00c57c1a94b226624bc",
    )

    # Compatible Python versions
    depends_on("python@3.7:", type=("build", "run"))

    # Build dependencies
    depends_on("py-setuptools", type="build")

    def config_settings(self, spec, prefix):
        settings = {}
        return settings
