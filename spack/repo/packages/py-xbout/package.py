# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class PyXbout(PythonPackage):
    """xBOUT provides an interface for collecting the output data from a BOUT++ simulation into an xarray dataset in an efficient and scalable way, as well as accessor methods for common BOUT++ analysis and plotting tasks."""

    homepage = "https://github.com/boutproject/xBOUT"
    pypi = "xbout/xbout-0.3.7.tar.gz"

    # Set a maintainer if submitting this package to the spack repo
    # maintainers("github_user1", "github_user2")

    license("Apache-2.0")

    version(
        "0.3.7",
        sha256="51b6bcc888553037a623f68dccfe7755ca409801d5b2dd1b8b1ecaca78c1eff1",
    )

    # Compatible Python versions
    depends_on("python@3.9:", type=("build", "run"))

    # Build dependencies
    depends_on("py-setuptools@65:", type="build")
    depends_on("py-setuptools-scm@7:+toml", type="build")
    depends_on("py-wheel@0.29.0:", type="build")

    # Runtime dependencies
    depends_on("py-animatplot-ng@0.4.2:", type=("build", "run"))
    depends_on("py-boutdata@0.1.4:", type=("build", "run"))
    depends_on("py-dask@2.10.0:+array", type=("build", "run"))
    depends_on("py-gelidum@0.5.3:", type=("build", "run"))
    depends_on("py-matplotlib@3.3.3:", type=("build", "run"))
    depends_on("py-natsort@5.5.0:", type=("build", "run"))
    depends_on("py-netcdf4@1.4.0:", type=("build", "run"))
    depends_on("py-pillow@6.1.0:", type=("build", "run"))
    depends_on("py-xarray@2023.01.0:", type=("build", "run"))

    def config_settings(self, spec, prefix):
        settings = {}
        return settings
