Hermes-3 manual
===============

This manual is available on [ReadTheDocs](https://hermes3.readthedocs.io/en/latest/).

The manual is written in
[ReStructuredText](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html).
To build it into a PDF or HTML output, install
[sphinx-build](https://www.sphinx-doc.org/en/master/usage/installation.html) e.g. with pip:
    pip install sphinx

Part of the documentation is auto-generated from the source code using Doxygen.
To enable Doxygen capability, install Doxygen:

    sudo apt install doxygen

And install Breathe, which connects Doxygen and Sphinx:

    pip install breathe

NOTE: Sphinx will happily compile the documentation without Doxygen as long as you don't have Breathe installed. 

To compile the documentation you first have to run Doxygen to generate the hermes-3/docs/doxygen/hermes3 directory:

    cd hermes-3/docs/doxygen
    Doxygen doxyfile

And finally compile the html docs with Sphinx:

    cd hermes-3/docs
    sphinx-build sphinx build

Where "sphinx" and "build" are the source and build directories, respectively.

