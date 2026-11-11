#!/bin/bash

set -e

if [[ ! -f $HOME/local_bout/.bout_install_complete ]]; then
    echo "****************************************"
    echo "Building BOUT++"
    echo "****************************************"
    if [[ -d $HOME/local_bout ]]; then
        echo ""
        echo "Cleaning up corrupted installation"
        echo ""
        rm -Rf "$HOME/local_bout"
    fi
    mkdir -p external/BOUT-dev/build/ && cd external/BOUT-dev/build/
    cmake -DCMAKE_INSTALL_PREFIX="$HOME/local_bout" \
          -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
          $BOUT_CONFIGURE_OPTIONS \
          ..
    make  -j4 && make -j4 install && touch "$HOME/local_bout/.bout_install_complete"
    echo "****************************************"
    echo "Finished building BOUT++"
    echo "****************************************"
else
    echo "****************************************"
    echo "BOUT++ already installed"
    echo "****************************************"
fi
