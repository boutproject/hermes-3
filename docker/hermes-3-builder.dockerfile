# Build as "hermes-3-builder"
# with sudo docker build -f docker/hermes-3-builder.dockerfile -t hermes-3-builder .

# Use a spack image with a pinned SHA
FROM spack/ubuntu-jammy@sha256:d9acf9ed998cbde8d12bd302c5921291086bfe6f70d2d0e26908fdc48c272324 AS builder
# Make sure that spack is available if we need to launch a terminal in the image
RUN spack_path=$(which spack) && \
    echo "alias spack=\"${spack_path}\"" >> /root/.bashrc

# Copy in the global spack configuration file
COPY docker/image_ingredients/spack_config.yaml /root/.spack/config.yaml

# Update package lists
RUN apt-get -yqq update && apt-get -yqq upgrade \
 && apt-get -yqq install --no-install-recommends ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Let spack use binary caches, which can significantly speed up the build process
RUN spack mirror add spack-cache https://cache.spack.io/
RUN spack mirror add develop https://binaries.spack.io/develop
RUN spack buildcache keys --install --trust

# Install GCC. You can set the version via build-args
ARG GCC_VERSION=14.2.0
# Building GCC from scratch takes a very long time. Instead, we
# force spack to use buildcache, which is much faster. If this
# fails, make sure that the requested GCC_VERSION is available
# in the build cache
RUN spack install -v --use-buildcache only gcc@${GCC_VERSION}
# Add this compiler to list of available compilers
RUN spack load gcc@${GCC_VERSION} && spack compiler find

# Make a new spack environment folder
RUN mkdir -p /opt/spack-environment
WORKDIR /opt/spack-environment
# Copy in the spack.yaml file
COPY docker/image_ingredients/spack.yaml /opt/spack-environment/spack.yaml
# Add the gcc version we just added to the environment and install
RUN spack env activate . && spack add gcc@${GCC_VERSION} && spack install -v --fail-fast

# # Make an 'entrypoint.sh' script which activates the spack environment
RUN spack env activate --sh -d . > activate.sh

# # Copy in a script which runs before any instance is launched
COPY docker/image_ingredients/docker_entrypoint.sh /entrypoint.sh
RUN chmod a+x /entrypoint.sh

# Set the default entrypoint and command for the image
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash"]