# Build as "hermes-3-builder"
# with sudo docker build -f docker/hermes-3-builder.dockerfile -t hermes-3-builder .

# Use a spack image pinned to the 1.2 release. The reference is tag@digest:
# Docker resolves *solely* by the @sha256 digest (immutable - the exact bytes
# never change even if the tag is re-pushed) and does NOT verify the tag against
# it. The :1.2.0 tag is unverified documentation only, so keep it in sync with
# the digest by hand when bumping.
# N.B. The spack 1.1 release has a bug in the PETSc package that causes this build to fail; 1.2 fixes it.
FROM spack/ubuntu-noble:1.2.0@sha256:dafccd1a2e77c61b5b6f81c06bbdabd999e6886daa405b3ac9011c5e4d98f8fa AS builder

# Location of a self-hosted OCI Spack build cache (see below) and the registry
# user to authenticate as. When SPACK_OCI_USER is empty (e.g. local builds or
# fork PRs with no token) the OCI cache is skipped and only the public mirror
# is used.
ARG SPACK_OCI_CACHE=oci://ghcr.io/boutproject/hermes-3-spack-cache
ARG SPACK_OCI_USER=""
ENV SPACK_OCI_CACHE=${SPACK_OCI_CACHE} \
    SPACK_OCI_USER=${SPACK_OCI_USER}

# Make sure that spack is available if we need to launch a terminal in the image
RUN spack_path=$(which spack) && \
    echo "alias spack=\"${spack_path}\"" >> /root/.bashrc

# Copy in the global spack configuration file
COPY docker/image_ingredients/spack_config.yaml /root/.spack/config.yaml

# No OS packages are installed here: the base image already ships the full GCC
# toolchain (gcc/g++/gfortran/make, libc6-dev) plus git, and cmake is built by
# Spack as part of the environment (listed in spack.yaml).
#
# Register gcc as an external package. It is already registered by the base
# image in /root/.spack/packages.yaml; re-running is harmless and acts as a
# safety net. In Spack 1.x compilers are graph nodes and the concretizer only
# accepts compilers found as external *packages* via `spack external find` (the
# legacy `spack compiler find`, which writes compilers.yaml, is ignored).
RUN spack external find gcc

# Check what architecture is supported by the builder
RUN /usr/lib64/ld-linux-x86-64.so.2 --help|grep supported

# Add Spack's public binary mirror so most dependencies can be downloaded as
# prebuilt binaries instead of compiled from source. Signed, so import the
# public signing keys and trust them. (Cache hits depend on the concretization
# matching what Spack built upstream, so this may only cover some packages.)
# In Spack 1.2 `buildcache keys --trust` prompts to confirm each key ([y/N]),
# which EOFs during a non-interactive build; --yes-to-all answers automatically.
RUN spack mirror add --scope site spack-public https://binaries.spack.io/develop \
 && spack buildcache keys --install --trust --yes-to-all

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir -p /opt/spack-environment
COPY docker/image_ingredients/spack.yaml /opt/spack-environment/spack.yaml

# Install the software. When registry credentials are provided (via the
# 'ghcr_token' build secret and SPACK_OCI_USER build-arg), also use a
# self-hosted OCI build cache: pull our own previously-compiled binaries
# before building, then push whatever we built back. Because the cache is
# keyed per-package hash, later builds only recompile packages that actually
# changed - even after edits to spack.yaml. Without a token, this is skipped
# and we fall back to the public mirror + compiling from source.
WORKDIR /opt/spack-environment
# Pin a hard, portable microarchitecture target for this arch. The
# `concretizer:targets:granularity: generic` clause in spack.yaml is only a
# *preference*: with reuse enabled (the default, and required for the OCI cache
# below) the concretizer happily reuses a prebuilt binary with a newer, native
# microarch - e.g. armv9.0a from the public mirror, or x86_64_v4 (AVX-512) built
# on a GitHub runner - which then SIGILLs ("Illegal instruction") on older
# hardware or under a Docker VM. A `require: target=...` is a hard constraint
# that reuse cannot override, so we inject it per-arch here (the shared
# spack.yaml cannot branch on architecture). aarch64/x86_64_v3 are the portable
# baselines; an external gcc can always target an older ISA than its build host,
# so this does not make the compiler node unsatisfiable.
ARG HERMES_TARGET=x86_64_v4
RUN spack -e . config add "packages:all:require:target=${HERMES_TARGET}"
# Save progress even on failure. Without --fail-fast, Spack builds every package
# whose dependencies succeeded before giving up, so a late failure still leaves
# most packages installed. We then push whatever installed (regardless of the
# install exit code) so the next build reuses it as a prebuilt binary, and
# finally re-propagate the exit code so a failed build still fails CI.
RUN --mount=type=secret,id=ghcr_token \
    spack env activate . && \
    if [ -n "${SPACK_OCI_USER}" ] && [ -s /run/secrets/ghcr_token ]; then \
      export SPACK_OCI_TOKEN="$(cat /run/secrets/ghcr_token)" && \
      spack mirror add --unsigned \
        --oci-username-variable SPACK_OCI_USER \
        --oci-password-variable SPACK_OCI_TOKEN \
        hermes-oci "${SPACK_OCI_CACHE}" && \
      { spack install ; rc=$? ; \
        spack buildcache push --unsigned --update-index --without-build-dependencies hermes-oci || true ; \
        exit $rc ; } ; \
    else \
      spack install --fail-fast ; \
    fi

# Make an 'entrypoint.sh' script which activates the spack environment
RUN spack env activate --sh -d . > activate.sh

# Copy in a script which runs before any instance is launched
COPY docker/image_ingredients/docker_entrypoint.sh /entrypoint.sh
RUN chmod a+x /entrypoint.sh

# Set the default entrypoint and command for the image
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash"]
