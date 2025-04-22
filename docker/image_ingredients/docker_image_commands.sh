#!/bin/bash

# Define color codes for pretty output
LIGHTRED='\033[1;31m'
LIGHTGREEN='\033[1;32m'
LIGHTBLUE='\033[1;34m'
NOCOLOR='\033[0m' # Reset to no color

# Define function for printing messages in different colors
print_message() {
    local color="$1"
    local message="$2"
    printf "\n${color}%s${NOCOLOR}\n" "$message"
}

# Simplified function calls for different message types
warn() { print_message "${LIGHTRED}" "$1"; }
notice() { print_message "${LIGHTGREEN}" "$1"; }
info() { print_message "${LIGHTBLUE}" "$1"; }
status() { print_message "" "$1"; } # No color for status messages

# First argument is interpreted as the mode
executable_name=$0
mode=$1
run_directory=$2
# Mode must be one of the handled options
valid_modes=("build_boutpp" "build_hermes" "build_both" "run" "build_hermes_and_run" "build_both_and_run" "fix_permissions")
if [[ " ${valid_modes[@]} " =~ " ${mode} " ]]; then
    if [[ -n "$run_directory" ]]; then
        notice "Running ${executable_name} in mode '${mode}' with run directory '${run_directory}'"
    else
        notice "Running ${executable_name} in mode '${mode}'"
    fi
else
    warn "Error: invalid mode '$mode' for Hermes build command."
    warn "Valid modes are: ${valid_modes[*]}" >&2
    exit 1 # Return a non-zero exit code to indicate failure
fi

build_boutpp () {
    local source_dir="${BOUTPP_SRC_DIR}"
    local config_file="${BOUTPP_CONFIG}"

    if [[ -d "${BOUTPP_SRC_DIR_OVERRIDE}" ]]; then
        source_dir="${BOUTPP_SRC_DIR_OVERRIDE}"
    fi
    if [[ -f "${BOUTPP_CONFIG_OVERRIDE}" ]]; then
        config_file="${BOUTPP_CONFIG_OVERRIDE}"
    fi
    info "Removing previous build directory $BOUTPP_BUILD_DIR"
    rm -rf ${BOUTPP_BUILD_DIR}
    notice "Configuring BOUT++ from ${source_dir} using the config file ${config_file}"
    cmake  -Wno-dev -B "${BOUTPP_BUILD_DIR}" -S "${source_dir}" -C "${config_file}"
    [ $? -ne 0 ] && warn "BOUT++ configuration failed" && exit 1
    notice "Finished configuring BOUT++. Starting build"
    cmake --build "${BOUTPP_BUILD_DIR}" --parallel
    [ $? -ne 0 ] && warn "BOUT++ build failed" && exit 1
    notice "Finished building BOUT++"
}

build_hermes () {
    local source_dir="${HERMES_SRC_DIR}"
    local config_file="${HERMES_CONFIG}"

    if [[ -d "${HERMES_SRC_DIR_OVERRIDE}" ]]; then
        source_dir="${HERMES_SRC_DIR_OVERRIDE}"
    fi
    if [[ -f "${HERMES_CONFIG_OVERRIDE}" ]]; then
        config_file="${HERMES_CONFIG_OVERRIDE}"
    fi
    info "Removing previous build directory $HERMES_BUILD_DIR"
    rm -rf ${HERMES_BUILD_DIR}
    notice "Configuring Hermes-3 from ${source_dir} using the config file ${config_file}"
    cmake  -Wno-dev -B "${HERMES_BUILD_DIR}" -S "${source_dir}" -C "${config_file}" \
           -DCMAKE_PREFIX_PATH=${BOUTPP_BUILD_DIR}
    [ $? -ne 0 ] && warn "Hermes-3 configuration failed" && exit 1
    notice "Finished configuring Hermes-3. Starting build"
    cmake --build "${HERMES_BUILD_DIR}" --parallel
    [ $? -ne 0 ] && warn "Hermes-3 build failed" && exit 1
    notice "Finished building Hermes-3"
}

run_hermes () {
    if [[ ! -n "$run_directory" ]]; then
        warn "Error: need to pass a directory to the 'build run' command."
        exit 1
    fi
    local run_dir="/hermes_project/${run_directory}"
    notice "Run directory on host machine: ${HOST_DIR}/${run_directory}"
    notice "Run directory in image: ${run_dir}"
    if [[ ! -d "${run_dir}" ]]; then
        warn "Error: run directory '${run_dir}' does not exist in the work directory on image."
        exit 1
    elif [[ ! -f "${run_dir}/BOUT.inp" ]]; then
        warn "Error: run directory '${run_dir}' does not contain a BOUT.inp file."
        exit 1
    else
        notice "Running Hermes-3"
        /hermes_project/build/hermes-3-build/hermes-3 -d $run_dir
    fi
}

fix_permissions () {
    notice "Adjusting permissions of /hermes_project/work for hermes_user"
    chown -R ${PUID}":"${PGID} /hermes_project/work
}

case "$mode" in
    "fix_permissions")
    fix_permissions
    ;;
    "build_boutpp")
    build_boutpp
    ;;
    "build_hermes")
    build_hermes
    ;;
    "build_both")
    build_boutpp && build_hermes
    ;;
    "run")
    run_hermes
    ;;
    "build_hermes_and_run")
    build_hermes && run_hermes
    ;;
    "build_both_and_run")
    build_boutpp && build_hermes && run_hermes
esac