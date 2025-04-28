from spack import *


class PyXhermes(PythonPackage):
    """xhermes Python package"""

    homepage = "https://github.com/boutproject/xhermes"
    url = "https://github.com/boutproject/xhermes"
    git = "https://github.com/boutproject/xhermes"

    # maintainers = ["Mike Kryjak"]

    version("main", branch="main", no_cache=True)

    # Compatible Python versions
    depends_on("python@3.6:", type=("build", "run"))

    # Build dependencies
    depends_on("py-setuptools", type="build")

    # Runtime dependencies
    depends_on("py-xbout", type=("build", "run"))
    depends_on("py-xarray", type=("build", "run"))
