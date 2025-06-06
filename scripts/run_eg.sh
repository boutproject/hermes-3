#!/bin/bash
#
# A convenience script to run examples.
# Run without any arguments for usage instructions.


REPO_ROOT=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" )/..)" &> /dev/null && pwd )
EGS_ROOT="${REPO_ROOT}/examples"
RUNS_ROOT="${REPO_ROOT}/runs"

#--------------------------------------------------------------------------------------------------
# Helper functions
echo_usage() {
    echo "Usage:"
    echo "    $0 [example_relpath] <-n num_MPI> <-b build_dir>"
    echo "    where"
    echo "       example_relpath is the path to the example (sub)directory relative to ./examples"
    echo "       num_MPI is the number of mpi processes to use"
    echo "       build_dir is the directory containing the hermes-3 exec (Defaults to one directory above this script)."
    echo
    echo "    e.g. $0 tokamak-1D/1D-threshold -n 8 -b builds/my_build"
}

# Execute a command, echoing it back first
execute() {
    local run_cmd=$1
    echo "----------------------------------------------------------------------------------------------------"
    echo "Executing [$run_cmd]"
    eval "$run_cmd"
    echo "----------------------------------------------------------------------------------------------------"
}

# Create a run directory for the example in $RUNS_ROOT
generate_run_dir() {
    local eg_dir="$1"
    local run_dir="$2"
    run_dir="${RUNS_ROOT}/$eg_relpath"
    if [ -e "$run_dir" ]; then
        read -p "Overwrite existing run directory at $run_dir? (Y/N): " choice && [[ $choice == [yY] || $choice == [yY][eE][sS] ]] || exit 5
        \rm -rf "$run_dir"
    fi
    mkdir -p "$(dirname $run_dir)"
    cp -r "$eg_dir" "$run_dir"
}

# Look for hermes-3 exec in the root of the repo by default
set_default_build_dir() {
    build_dir="$REPO_ROOT/"
}

# Parse CL args
parse_args() {
    if [ $# -lt 1 ]; then
        echo_usage
        exit 1
    fi
    POSITIONAL_ARGS=()
    while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--build-dir)
        build_dir=$(realpath "$2")
        shift 2
        ;;
        -n|--num_mpi)
        nmpi="$2"
        shift 2
        ;;
        -*)
        echo "Unknown option $1"
        exit 2
        ;;
        *)
        # Save positional args in an array
        POSITIONAL_ARGS+=("$1")
        shift
        ;;
    esac
    done

    # Restore and extract positional args
    set -- "${POSITIONAL_ARGS[@]}"
    eg_relpath=$1
}

# echo back options
report_options() {
    echo "Options:"
    echo "          example relpath : $eg_relpath"
    echo "                    n MPI : $nmpi"
    echo ""
}

# Check that the executable and example directory both exist
validate_paths() {
    local exec=$1
    local eg_dir=$2
    if [ ! -f "$exec" ]; then
        echo "No executable found at $exec"
        exit 3
    fi
    if [ ! -d "$eg_dir" ]; then
        echo "No example directory found at $eg_dir"
        exit 4
    fi
}
#--------------------------------------------------------------------------------------------------

# Default options
eg_relpath='Not set'
nmpi='4'
build_dir='Not set'
set_default_build_dir

# Parse command line args and report resulting options
parse_args $*
report_options

# Set paths to the executable and example directory
h3_exec="$build_dir/hermes-3"
eg_dir="$EGS_ROOT/$eg_relpath"

# Validate exec, examples paths
validate_paths "$h3_exec" "$eg_dir"

# Set up a run directory, confirming overwrite if it already exists
run_dir="${RUNS_ROOT}/$eg_relpath"
generate_run_dir "$eg_dir" "$run_dir"

# Set up run_cmd and execute in run_dir
run_cmd="mpirun -n $nmpi $h3_exec -d $run_dir"
execute "$run_cmd"
