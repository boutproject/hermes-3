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

# Helper function to check exit codes and exit on failure
check_exit_code() {
    local exit_code=$?
    local message=$1
    if [[ $exit_code -ne 0 ]]; then
        warn "${message}"
        exit $exit_code
    fi
}

# First argument is interpreted as the mode
executable_name=$0
mode=$1
run_directory=$2
# Optional third argument: number of MPI ranks for the "run" modes (default 1)
nprocs=$3
# Mode must be one of the handled options
valid_modes=("build_boutpp" "build_hermes" "build_both" "run" "build_hermes_and_run" "build_both_and_run" "fix_permissions")

# Validate mode
mode_is_valid=false
for valid_mode in "${valid_modes[@]}"; do
    if [[ "$mode" == "$valid_mode" ]]; then
        mode_is_valid=true
        break
    fi
done

if ! $mode_is_valid; then
    warn "Error: invalid mode '$mode' for Hermes build command."
    warn "Valid modes are: ${valid_modes[*]}" >&2
    exit 1 # Return a non-zero exit code to indicate failure
fi

# Refuse to run as root unless the mode genuinely needs it. Every mode except
# fix_permissions writes build/run output into the mounted work/ volume; doing
# that as root leaves root-owned files on the host (compose runs these services
# as user: ${UID}:${GID} precisely to avoid that). fix_permissions is exempt
# because its whole job is to chown that volume. Set HERMES_ALLOW_ROOT=1 to
# override, e.g. when the host user genuinely is root.
root_allowed_modes=("fix_permissions")
mode_allows_root=false
for allowed_mode in "${root_allowed_modes[@]}"; do
    if [[ "$mode" == "$allowed_mode" ]]; then
        mode_allows_root=true
        break
    fi
done

if [[ "$(id -u)" -eq 0 ]] && ! $mode_allows_root && [[ "${HERMES_ALLOW_ROOT}" != "1" ]]; then
    warn "Error: refusing to run mode '$mode' as root."
    warn "It would create root-owned files in the mounted work/ directory on the host."
    warn "Run as your host user (compose sets user: \${UID}:\${GID}), or set HERMES_ALLOW_ROOT=1 to override."
    exit 1
fi

# Announce mode
if [[ -n "$run_directory" ]]; then
    notice "Running ${executable_name} in mode '${mode}' with run directory '${run_directory}'"
else
    notice "Running ${executable_name} in mode '${mode}'"
fi

build_boutpp () {
    local source_dir="${BOUTPP_SRC_DIR}"
    local config_file="${BOUTPP_CONFIG}"
    local build_dir="${BOUTPP_BUILD_DIR_OVERRIDE}"

    # Check for overrides
    [[ -d "${BOUTPP_SRC_DIR_OVERRIDE}" ]] && source_dir="${BOUTPP_SRC_DIR_OVERRIDE}"
    [[ -f "${BOUTPP_CONFIG_OVERRIDE}" ]] && config_file="${BOUTPP_CONFIG_OVERRIDE}"

    if [[ -d "${build_dir}" ]]; then
        info "Removing previous build directory ${build_dir}"
        rm -rf "${build_dir}"
    else
        info "Build directory ${build_dir} is empty."
    fi

    notice "Configuring BOUT++ from ${source_dir} using the config file ${config_file}"
    cmake -Wno-dev -B "${build_dir}" -S "${source_dir}" -C "${config_file}"
    check_exit_code "BOUT++ configuration failed"

    local jobs="${HERMES_BUILD_JOBS:-4}"
    notice "Finished configuring BOUT++. Starting build with --parallel ${jobs} (set via HERMES_BUILD_JOBS)"
    cmake --build "${build_dir}" --parallel "${jobs}"
    check_exit_code "BOUT++ build failed"

    notice "Finished building BOUT++"
}

build_hermes () {
    local source_dir="${HERMES_SRC_DIR}"
    local config_file="${HERMES_CONFIG}"
    local boutpp_dir="${BOUTPP_BUILD_DIR}"
    local build_dir="${HERMES_BUILD_DIR_OVERRIDE}"

    # Check for overrides
    [[ -d "${HERMES_SRC_DIR_OVERRIDE}" ]] && source_dir="${HERMES_SRC_DIR_OVERRIDE}"
    [[ -f "${HERMES_CONFIG_OVERRIDE}" ]] && config_file="${HERMES_CONFIG_OVERRIDE}"
    [[ -d "${BOUTPP_BUILD_DIR_OVERRIDE}" ]] && boutpp_dir="${BOUTPP_BUILD_DIR_OVERRIDE}"

    if [[ -d "${build_dir}" ]]; then
        info "Removing previous build directory ${build_dir}"
        rm -rf "${build_dir}"
    else
        info "Build directory ${build_dir} is empty."
    fi

    notice "Configuring Hermes-3 from ${source_dir} using the config file ${config_file}"
    cmake -Wno-dev -B "${build_dir}" -S "${source_dir}" -C "${config_file}" \
           -DCMAKE_PREFIX_PATH="${boutpp_dir}"
    check_exit_code "Hermes-3 configuration failed"

    local jobs="${HERMES_BUILD_JOBS:-4}"
    notice "Finished configuring Hermes-3. Starting build with --parallel ${jobs} (set via HERMES_BUILD_JOBS)"
    cmake --build "${build_dir}" --parallel "${jobs}"
    check_exit_code "Hermes-3 build failed"

    notice "Finished building Hermes-3"
}

run_hermes () {
    if [[ -z "$run_directory" ]]; then
        warn "Error: need to pass a non-empty directory to the run command."
        exit 1
    fi

    # Number of MPI ranks. Defaults to a single process; must be a positive integer.
    local ranks="${nprocs:-1}"
    if ! [[ "$ranks" =~ ^[1-9][0-9]*$ ]]; then
        warn "Error: number of processes must be a positive integer, got '${ranks}'."
        exit 1
    fi

    local run_dir="/hermes_project/${run_directory}"
    local hermes_dir="${HERMES_BUILD_DIR}"

    [[ -d "${HERMES_BUILD_DIR_OVERRIDE}" ]] && hermes_dir="${HERMES_BUILD_DIR_OVERRIDE}"
    info "Using ${hermes_dir}/hermes-3"

    notice "Run directory on host machine: ${HOST_DIR}/${run_directory}"
    notice "Run directory in image: ${run_dir}"
    if [[ ! -d "${run_dir}" ]]; then
        warn "Error: run directory '${run_dir}' does not exist in the work directory on image."
        exit 1
    elif [[ ! -f "${run_dir}/BOUT.inp" ]]; then
        warn "Error: run directory '${run_dir}' does not contain a BOUT.inp file."
        exit 1
    elif [[ "$ranks" -eq 1 ]]; then
        notice "Running Hermes-3 in ${run_dir} on a single process"
        # Execute hermes-3 directly without changing directory
        "${hermes_dir}/hermes-3" -d "${run_dir}"
        check_exit_code "Hermes-3 execution failed"
        notice "Finished running Hermes-3"
    else
        # Only pass --oversubscribe when the request exceeds the cores Docker
        # exposes. Within the core count, mpirun runs without it (so a genuine
        # over-request elsewhere still fails loudly); above it we warn, since
        # oversubscribed ranks time-slice the same CPU and run slowly.
        local ncores
        ncores=$(nproc 2>/dev/null || echo 1)
        local oversubscribe=""
        if [[ "$ranks" -gt "$ncores" ]]; then
            warn "Warning: requested ${ranks} ranks but only ${ncores} core(s) are available to the container."
            warn "Running with --oversubscribe; this will be slow. Raise Docker's CPU allocation for better performance."
            oversubscribe="--oversubscribe"
        fi
        notice "Running Hermes-3 in ${run_dir} on ${ranks} processes with mpirun"
        mpirun ${oversubscribe} -np "${ranks}" "${hermes_dir}/hermes-3" -d "${run_dir}"
        check_exit_code "Hermes-3 execution failed"
        notice "Finished running Hermes-3"
    fi
}

fix_permissions () {
    notice "Adjusting permissions of /hermes_project/work for hermes_user (${PUID}:${PGID})"
    chown -R "${PUID}:${PGID}" /hermes_project/work
    check_exit_code "Failed to adjust permissions on /hermes_project/work"
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
