#!/bin/bash

set -e

export PATH=${HOME}/.local/bin:${PATH}
pip3 install --user --upgrade pip setuptools
pip3 install --user --upgrade scipy numpy
for package in $@
do
    if test $package == "cython"
    then
        # fast install Cython
        pip3 install --user Cython --install-option="--no-cython-compile"
    elif test $package == "something_else"
    then
        pip3 install what_we_need
    else
        pip3 install --user $package
    fi
done

# Install Hermes-3
pip install --user git+https://github.com/boutproject/xhermes
