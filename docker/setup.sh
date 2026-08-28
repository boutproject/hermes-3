#!/bin/bash
mkdir -p work

echo "UID=$(id -u)" > .env
echo "GID=$(id -g)" >> .env
echo "HOST_DIR=$PWD" >> .env
echo "HERMES_TAG=latest" >> .env
echo "JUPYTER_TAG=latest" >> .env
echo "HERMES_BUILD_JOBS=4" >> .env
