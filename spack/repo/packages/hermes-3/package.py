from spack import *


class Hermes3(CMakePackage):
    """Hermes3"""

    homepage = "https://github.com/bendudson/hermes-3"

    git = "https://github.com/bendudson/hermes-3.git"
    version("working", branch="master")
    version("master", branch="master")

    variant("petsc", default=False, description="Builds with PETSc support.")
    variant("sundials", default=True, description="Builds with SUNDIALS support.")
    variant("xhermes", default=True, description="Builds xhermes.")

    depends_on("cmake@3.24:", type="build")
    # depends_on("adios2", type=("build", "link", "run"))
    depends_on("fftw", type=("build", "link", "run"))
    depends_on("mpi", type=("build", "link", "run"))
    depends_on("netcdf-cxx4", type=("build", "link", "run"))
    depends_on("py-cython", type=("build", "link", "run"))
    depends_on("py-jinja2", type=("build", "link", "run"))
    depends_on("py-netcdf4", type=("build", "link", "run"))
    depends_on("py-xhermes", type=("build", "link", "run"))

    # Dependencies controlled by variants
    depends_on(
        "petsc+hypre+mpi~debug~fortran", when="+petsc", type=("build", "link", "run")
    )
    # Could set up Sundials as a spack dependency here; just downloaded it via the BOUT cmake flag for now
    # depends_on("sundials", when="+sundials", type=("build", "link", "run"))
    depends_on("py-xhermes", when="+xhermes", type=("build", "link", "run"))

    def cmake_args(self):
        args = []
        if "+petsc" in self.spec:
            args.append("-DBOUT_USE_PETSC=ON")
        if "+sundials" in self.spec:
            args.append("-DBOUT_DOWNLOAD_SUNDIALS=ON")
        return args
