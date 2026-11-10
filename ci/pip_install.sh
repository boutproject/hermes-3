#!/bin/bash

set -e

export PATH=${HOME}/.local/bin:${PATH}
# Upgrade to trusted version of pip, see https://discuss.python.org/t/information-for-python-users-regarding-emudhra-emsign-ca-certificate-security-incident/109576/4
pip3 install --user "pip @ https://files.pythonhosted.org/packages/3f/06/3e4420302363fbb71be3c7f35ce042f54bfca330e3e66b3a2e3c4cca579f/pip-26.2.2-py3-none-any.whl#sha256=fd7f4dfa8190e014013b9f08f3e973f4e4c26ffd179a204895533a7254813a70"
# From there can upgrade to latest pip
pip3 install --user --upgrade pip setuptools
pip3 install --user -r ci/requirements.txt
