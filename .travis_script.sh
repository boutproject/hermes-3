#!/bin/bash

# Exit on error
set -e

RED_FG="\033[031m"
RESET_FG="\033[039m"

echo "****************************************"
echo "Fetching BOUT++ master"
echo "****************************************"

git clone --depth 1 https://github.com/boutproject/BOUT-dev.git
cd BOUT-dev

export MAKEFLAGS="-j 2 -k"
echo "****************************************"
echo "Configuring with $BOUT_CONFIGURE_OPTIONS"
echo "****************************************"
conf=0
time ./configure $BOUT_CONFIGURE_OPTIONS MAKEFLAGS="$MAKEFLAGS" || conf=$?
if test $conf -gt 0
then
    echo -e $RED_FG
    echo "**************************************************"
    echo "Printing config.log:"
    echo "**************************************************"
    echo -e $RESET_FG
    echo
    cat config.log
    echo
    echo -e $RED_FG
    echo "**************************************************"
    echo "Printing config-build.log:"
    echo "**************************************************"
    echo -e $RESET_FG
    echo
    cat config-build.log
    exit $conf
fi
export PYTHONPATH=$(pwd)/tools/pylib/:$PYTHONPATH

for target in ${MAIN_TARGET[@]}
do
    make_exit=0
    time make $target || make_exit=$?
    if [[ $make_exit -gt 0 ]]; then
	make clean > /dev/null
	echo -e $RED_FG
	echo "**************************************************"
	echo "Printing make commands:"
	echo "**************************************************"
	echo -e $RESET_FG
	echo
	make -n $target
	exit $make_exit
    fi
done

# Change to parent directory (Hermes-3)
cd ..

echo -e $RED_FG
echo "**************************************************"
echo "Compiling Hermes-3"
echo "**************************************************"
echo -e $RESET_FG

time make BOUT_TOP=BOUT-dev

echo -e $RED_FG
echo "**************************************************"
echo "Running tests"
echo "**************************************************"
echo -e $RESET_FG

# Run tests
time make check

