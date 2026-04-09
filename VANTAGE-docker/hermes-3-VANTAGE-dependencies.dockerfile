# Build stage with Spack pre-installed and ready to be used
FROM spack/ubuntu-noble@sha256:a1c7d1dcfea874f74ec827851f04fc116ca54869c30a766c92c5299967a4f33c
# above FROM equivalent to
# FROM spack/ubuntu-noble:1.1.0
# the image sha256 hash is obtained from the following url
# https://hub.docker.com/layers/spack/ubuntu-noble/1.1.0/images/sha256-a1c7d1dcfea874f74ec827851f04fc116ca54869c30a766c92c5299967a4f33c
# we use the sha256 hash to pin to a specific version
# in case the image associated with tag :1.1.0 is updated

RUN apt update && \
 apt install -y git --no-install-recommends \
 && rm -rf /var/lib/apt/lists/*
# use develop version of spack repos
RUN sed -i '/^[[:space:]]*branch:/ s|releases/v2025\.11|develop|g' /opt/spack/etc/spack/defaults/base/repos.yaml
# update spack repos
RUN spack repo update
# find the gcc compiler
RUN spack compiler find
# copy the files in the hermes-3 repo needed to install dependencies
COPY ./external /root/hermes-3/external
COPY ./spack.yaml /root/hermes-3/spack.yaml
# set the workdir
WORKDIR /root/hermes-3

# Install hermes-3 dependencies via spack
# Activate the hermes-3 environment, install dependencies
RUN <<EOF
# use generic microarchitectures to avoid conflicts between runners using zen2, icelake, etc
sed -i 's/^[[:space:]]*granularity:[[:space:]]*microarchitectures/      granularity: generic/' spack.yaml
# demand target=x86_64_v3, assume that spack.yaml has unique block
# all:
#   providers:
#     mpi: [mpich, openmpi ]
# which we mutate to
# all:
#   providers:
#     require: "target=x86_64_v3"
#     mpi: [mpich, openmpi ]
sed -i '/^[[:space:]]*providers:/i\      require: "target=x86_64_v3"' spack.yaml
spack env activate . -v gcc
# Install the dependencies
spack install -j 4 --only dependencies
# Uninstall the top-level packages that we expect to develop regularly
spack uninstall -y --dependents boutpp@develop || true
spack uninstall -y --dependents vantagereactions@working || true
spack uninstall -y --dependents neso-particles@working || true
EOF
# set workdir back to /root
WORKDIR /root
# remove hermes-3 clone, leaving dependencies installed
RUN rm -rf /root/hermes-3

ENTRYPOINT ["/bin/bash"]
