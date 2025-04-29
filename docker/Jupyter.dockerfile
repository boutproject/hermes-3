# Build this as "hermes3-jupyter"
FROM jupyter/scipy-notebook

RUN git clone https://github.com/boutproject/xhermes && cd xhermes && pip install -e .