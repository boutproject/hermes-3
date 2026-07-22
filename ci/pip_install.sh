#!/bin/bash

set -e

export PATH=${HOME}/.local/bin:${PATH}
pip3 install --user --upgrade pip setuptools
pip3 install --user -r ci/requirements.txt
