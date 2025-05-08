# Build as "hermes-3-builder"
# with sudo docker build -f docker/hermes-3-builder.dockerfile -t hermes-3-builder .

# Use a spack image with a pinned SHA
FROM spack/ubuntu-jammy:develop AS builder
# Make sure that spack is available if we need to launch a terminal in the image
RUN spack_path=$(which spack) && \
    echo "alias spack=\"${spack_path}\"" >> /root/.bashrc

# Tell spack to install into /opt/software
RUN spack config add "config:install_tree:root:'/opt/software'"

# Install GCC. You can set the version via build-args
ARG GCC_VERSION=14.2.0
RUN spack install gcc@${GCC_VERSION} && spack compiler add "$(spack location -i gcc@${GCC_VERSION})"

# Install the software. Add any dependencies you need here.
RUN mkdir -p /opt/spack-environment
WORKDIR /opt/spack-environment
RUN echo "Configuring" \
 && spack env create -d . \
 && spack env activate -d . \
 && spack config add "concretizer:unify:true" \
 && spack config add "view:'/opt/views/view'" \
 && spack config add "packages:all:prefer:[ '%gcc@${GCC_VERSION}', '^openmpi' ]" \
 && spack add gcc@${GCC_VERSION} \
              cmake \
              python@3.12 \
              py-cython \
              py-numpy \
              py-pip \
              py-jinja2 \
              netcdf-cxx4 \
              sundials \
              petsc+double+fftw+hdf5+hypre+mpi \
              slepc \
 && echo "Finished configuring. Concretizing." \
 && spack concretize \
 && echo "Finished concretizing. Installing." \
 && spack install --fail-fast

# Make an 'entrypoint.sh' script which activates the spack environment
RUN spack env activate --sh -d . > activate.sh
ENTRYPOINT [ "/bin/bash" ]
