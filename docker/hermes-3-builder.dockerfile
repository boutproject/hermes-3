# Build as "hermes-3-builder"
# with sudo docker build -f docker/hermes-3-builder.dockerfile -t hermes-3-builder .

# Spack 1.2 image, pinned by digest (the :1.2.0 tag is documentation only; keep
# it in sync by hand when bumping).
FROM spack/ubuntu-noble:1.2.0@sha256:dafccd1a2e77c61b5b6f81c06bbdabd999e6886daa405b3ac9011c5e4d98f8fa AS builder

# Self-hosted OCI Spack build cache and the registry user. When SPACK_OCI_USER
# is empty (local builds or fork PRs) the OCI cache is skipped and the stack
# builds from source.
ARG SPACK_OCI_CACHE=oci://ghcr.io/boutproject/hermes-3-spack-cache
ARG SPACK_OCI_USER=""
ENV SPACK_OCI_CACHE=${SPACK_OCI_CACHE} \
    SPACK_OCI_USER=${SPACK_OCI_USER}

# Extra flags for the no-credentials (source-build) `spack install`; empty in
# CI, local builds pass e.g. "--fresh --verbose".
ARG SPACK_INSTALL_EXTRA_ARGS=""

# GNU mirror(s) to seed readline patches from (space-separated, tried in order).
# ftpmirror.gnu.org, which Spack fetches these patches from, is unreachable from
# our build network. Set to "" to disable.
ARG GNU_MIRROR_URL="https://ftp.gnu.org/gnu https://mirrors.kernel.org/gnu"

# Make sure that spack is available if we need to launch a terminal in the image
RUN spack_path=$(which spack) && \
    echo "alias spack=\"${spack_path}\"" >> /root/.bashrc

# Copy in the global spack configuration file
COPY docker/image_ingredients/spack_config.yaml /root/.spack/config.yaml

# The base image already ships the GCC toolchain, git, and cmake-via-Spack, so
# no apt step is needed. Register gcc as an external package: in Spack 1.x
# compilers are graph nodes found only via `spack external find` (not the legacy
# `spack compiler find`).
RUN spack external find gcc

# The public Spack binary mirror is intentionally not added: some of its
# packages (e.g. python) are built with intel-oneapi, and reusing those leaks
# the Intel-only `-fp-model=strict` flag into gcc source builds. The stack
# builds from source with gcc instead.

RUN mkdir -p /opt/spack-environment
COPY docker/image_ingredients/spack.yaml /opt/spack-environment/spack.yaml

WORKDIR /opt/spack-environment
# Pin a hard, portable microarch target per-arch. spack.yaml's
# `granularity: generic` is only a preference: with reuse enabled (required by
# the OCI cache) the concretizer will happily reuse a newer native binary (e.g.
# armv9.0a or x86_64_v4/AVX-512) that then crashes with illegal instructions on
# older hardware. A `require: target=...` is a hard constraint reuse cannot
# override; we enforce it here because the shared spack.yaml can't branch on
# architecture.
# HERMES_TARGET overrides the microarch (CI builds compat variants, e.g.
# x86_64_v2); empty => portable default from `uname -m`.
ARG HERMES_TARGET=""
RUN if [ -z "${HERMES_TARGET}" ]; then \
      case "$(uname -m)" in \
        aarch64) HERMES_TARGET=aarch64 ;; \
        x86_64)  HERMES_TARGET=x86_64_v3 ;; \
        *)       HERMES_TARGET="$(uname -m)" ;; \
      esac ; \
    fi && \
    echo "Pinning Spack microarch target: ${HERMES_TARGET}" && \
    spack -e . config add "packages:all:require:target=${HERMES_TARGET}"

# Seed readline patches into a local Spack source mirror (see GNU_MIRROR_URL).
COPY docker/image_ingredients/seed_gnu_mirror.sh /tmp/seed_gnu_mirror.sh
RUN if [ -n "${GNU_MIRROR_URL}" ]; then \
      sh /tmp/seed_gnu_mirror.sh "${GNU_MIRROR_URL}" /opt/gnu-mirror && \
      spack mirror add --scope site local-gnu file:///opt/gnu-mirror ; \
    fi

# Install the environment. With OCI credentials, add the self-hosted cache,
# then install without --fail-fast and push whatever installed regardless of
# exit code (so a late failure still caches most packages for the next run),
# re-propagating the exit code so a failed build still fails CI. Without
# credentials, build from source (plus any SPACK_INSTALL_EXTRA_ARGS).
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
      spack install ${SPACK_INSTALL_EXTRA_ARGS} --fail-fast ; \
    fi

# Make an 'entrypoint.sh' script which activates the spack environment
RUN spack env activate --sh -d . > activate.sh

# Copy in a script which runs before any instance is launched
COPY docker/image_ingredients/docker_entrypoint.sh /entrypoint.sh
RUN chmod a+x /entrypoint.sh

# Set the default entrypoint and command for the image
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash"]
