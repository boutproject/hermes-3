# Build as "hermes-3-jupyter"
# with sudo docker build -f docker/hermes-3-jupyter.dockerfile -t hermes-3-jupyter .

# Maintained images now live at quay.io (docker.io/jupyter/scipy-notebook is frozen at Oct 2023).
# The reference is tag@digest: Docker resolves *solely* by the @sha256 digest
# (immutable - the exact bytes never change even if the tag is re-pushed) and does
# NOT verify the tag against it. The :2026-06-29 tag is unverified documentation
# only (jupyter images are date-tagged, not semver), so keep it in sync by hand.
FROM quay.io/jupyter/scipy-notebook:2026-06-29@sha256:52a7d9ee3faa90118d89db7729d6bb45db7cb030d9e2fdb096eb9979264d1ab6

RUN git clone https://github.com/boutproject/xhermes /home/jovyan/xhermes && cd /home/jovyan/xhermes && pip install -e .
USER root
RUN apt-get -yqq update && apt-get -yqq upgrade \
    && apt-get -yqq install --no-install-recommends doxygen \
    && rm -rf /var/lib/apt/lists/*
USER jovyan
RUN pip install pint-xarray sphinx sphinx_book_theme breathe

# Prevent an issue where only one process can open a HDF file at a time
ENV HDF5_USE_FILE_LOCKING=False
