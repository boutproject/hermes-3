# Build as "hermes-3-builder"
# with sudo docker build -f docker/hermes-3-builder.dockerfile -t hermes-3-builder .

# Use a spack image with a pinned SHA
FROM spack/ubuntu-jammy@sha256:d9acf9ed998cbde8d12bd302c5921291086bfe6f70d2d0e26908fdc48c272324 AS builder
# Make sure that spack is available if we need to launch a terminal in the image
RUN spack_path=$(which spack) && \
    echo "alias spack=\"${spack_path}\"" >> /root/.bashrc

# Copy in the global spack configuration file
COPY docker/image_ingredients/spack_config.yaml /root/.spack/config.yaml

# Install OS packages needed to build the software
RUN apt-get -yqq update && apt-get -yqq upgrade \
 && apt-get -yqq install --no-install-recommends git build-essential cmake \
 && rm -rf /var/lib/apt/lists/*
RUN spack compiler find && spack external find cmake

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