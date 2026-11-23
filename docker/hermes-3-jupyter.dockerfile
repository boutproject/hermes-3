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
# Make sure using a secure version of pip
RUN pip install "pip @ https://files.pythonhosted.org/packages/3f/06/3e4420302363fbb71be3c7f35ce042f54bfca330e3e66b3a2e3c4cca579f/pip-26.2.2-py3-none-any.whl#sha256=fd7f4dfa8190e014013b9f08f3e973f4e4c26ffd179a204895533a7254813a70"
RUN pip install --upgrade pip
RUN pip install pint-xarray sphinx sphinx_book_theme breathe

# Prevent an issue where only one process can open a HDF file at a time
ENV HDF5_USE_FILE_LOCKING=False
