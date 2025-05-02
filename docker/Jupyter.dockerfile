# Build as "hermes3-jupyter"
# with sudo docker build -f docker/Jupyter.dockerfile -t hermes3-jupyter docker

FROM jupyter/scipy-notebook

RUN git clone https://github.com/boutproject/xhermes /home/jovyan/xhermes && cd /home/jovyan/xhermes && pip install -e .
USER root
RUN apt-get -yqq update && apt-get -yqq upgrade \
    && apt-get -yqq install --no-install-recommends doxygen \
    && rm -rf /var/lib/apt/lists/*
USER jovyan
RUN pip install pint-xarray sphinx sphinx_book_theme breathe

# Prevent an issue where only one process can open a HDF file at a time
ENV HDF5_USE_FILE_LOCKING=False
