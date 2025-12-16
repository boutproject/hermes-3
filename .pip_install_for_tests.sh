#!/bin/bash

set -e

export PATH=${HOME}/.local/bin:${PATH}
pip3 install --user --upgrade pip~=20.0 setuptools~=46.1
pip3 install --user --upgrade scipy~=1.4 numpy~=1.18
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

# Install xBOUT to make sure we get latest master version
pip3 install  --user git+https://github.com/boutproject/xBOUT.git
# Install xHermes for Hermes-3 python
pip3 install  --user git+https://github.com/boutproject/xhermes.git

