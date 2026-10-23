#!/bin/bash

set -e

if [[ ! -d $HOME/local_bout/include/bout ]]; then
    echo "****************************************"
    echo "Building BOUT++"
    echo "****************************************"
    mkdir -p external/BOUT-dev/build/ && cd external/BOUT-dev/build/
    cmake -DCMAKE_INSTALL_PREFIX="$HOME/local_bout" \
          -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
          $BOUT_CONFIGURE_OPTIONS \
          ..
    make  -j4 && make -j4 install
    echo "****************************************"
    echo "Finished building BOUT++"
    echo "****************************************"
else
    echo "****************************************"
    echo "BOUT++ already installed"
    echo "****************************************"
fi
