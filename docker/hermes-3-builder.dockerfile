# Build as "hermes-3-builder"
# with sudo docker build -f docker/hermes-3-builder.dockerfile -t hermes-3-builder .

# Use a spack image with a pinned SHA - currently points to develop between the 1.1 and 1.2 releases.
# N.B. The spack 1.1 release has a bug in the PETSc package that causes this build to fail 
FROM spack/ubuntu-noble@sha256:d7784a53424fda1c528d8afe837841a6947e46b55fd4380779656d4b276f63a0 AS builder
# Make sure that spack is available if we need to launch a terminal in the image
RUN spack_path=$(which spack) && \
    echo "alias spack=\"${spack_path}\"" >> /root/.bashrc

# Copy in the global spack configuration file
COPY docker/image_ingredients/spack_config.yaml /root/.spack/config.yaml

# Install OS packages needed to build the software
RUN apt-get -yqq update && apt-get -yqq upgrade \
 && apt-get -yqq install --no-install-recommends git build-essential gfortran cmake \
 && rm -rf /var/lib/apt/lists/*
# Register gcc (c/cxx/fortran) and cmake as external packages. In Spack 1.x
# compilers are graph nodes and the concretizer only accepts compilers found
# as external *packages* via `spack external find`; the legacy `spack compiler
# find` (which writes compilers.yaml) is ignored, giving "no compilers found".
RUN spack external find gcc cmake

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir -p /opt/spack-environment
COPY docker/image_ingredients/spack.yaml /opt/spack-environment/spack.yaml

# Install the software
WORKDIR /opt/spack-environment
RUN spack env activate . && spack install --fail-fast

# Make an 'entrypoint.sh' script which activates the spack environment
RUN spack env activate --sh -d . > activate.sh

# Copy in a script which runs before any instance is launched
COPY docker/image_ingredients/docker_entrypoint.sh /entrypoint.sh
RUN chmod a+x /entrypoint.sh

# Set the default entrypoint and command for the image
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash"]