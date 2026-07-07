#!/bin/bash
mkdir -p work

echo "UID=$(id -u)" > .env
echo "GID=$(id -g)" >> .env
echo "HOST_DIR=$PWD" >> .env
echo "HERMES_TAG=latest" >> .env
echo "JUPYTER_TAG=latest" >> .env
echo "HERMES_BUILD_JOBS=4" >> .env

# Extra compiler flags appended when rebuilding in-container (build_hermes /
# build_boutpp). Leave empty to match the pulled image, which is already built
# for a fixed -march (a later -march here wins, so you can also retarget):
#   hermes-3:latest            -> -march=x86-64-v3  (armv8-a on arm64)
#   hermes-3:latest-x86_64_v2  -> -march=x86-64-v2
#   hermes-3:latest-x86_64_v1  -> -march=x86-64
# e.g. CMAKE_CXX_FLAGS="-march=native" or "-funroll-loops".
cat >> .env <<'EOF'
# See setup.sh for the image <-> -march correspondence.
CMAKE_CXX_FLAGS=
CMAKE_C_FLAGS=
EOF
