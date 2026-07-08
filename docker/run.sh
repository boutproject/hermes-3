#!/bin/bash
# Helper for running the Hermes-3 docker compose services.
#
# Usage:
#   Source this file and call run_docker, or run it directly:
#     ./run.sh <service> [arguments]
#
# Examples:
#   ./run.sh shell                        # interactive shell in the image
#   ./run.sh build_both                   # rebuild hermes and BOUT++
#   ./run.sh hermes work/test             # run a hermes-3 case (single rank)
#   ./run.sh hermes_parallel 4 work/test  # run a case on 4 MPI ranks
#   ./run.sh jupyter                      # start the Jupyter server
#   ./run.sh cleanup                      # tidy up orphaned/stopped containers
#   ./run.sh rm_docker                    # remove all images, containers and volumes

# Define color codes for pretty output
LIGHTRED='\033[1;31m'
LIGHTGREEN='\033[1;32m'
NOCOLOR='\033[0m' # Reset to no color

print_message() {
    local color="$1"
    local message="$2"
    printf "${color}%s${NOCOLOR}\n" "$message"
}

warn()   { print_message "${LIGHTRED}" "$1"; }
notice() { print_message "${LIGHTGREEN}" "$1"; }

# Services defined in docker-compose.yaml
VALID_SERVICES="shell sudo build_hermes build_boutpp build_both hermes fix_permissions jupyter"
# Services that require at least one argument (e.g. a case path under work/)
ARG_REQUIRED_SERVICES="hermes"

# Return 0 if $1 appears as a whitespace-separated word in $2
contains_word() {
  case " $2 " in
    *" $1 "*) return 0 ;;
    *) return 1 ;;
  esac
}

run_docker_help() {
  notice "Usage: run_docker <service> [arguments]"
  echo
  echo "Available services (from docker-compose.yaml):"
  echo "  shell            Interactive shell in the image"
  echo "  sudo             Interactive shell with root access"
  echo "  build_hermes     Rebuild hermes, using ./work/hermes-3 if available"
  echo "  build_boutpp     Rebuild BOUT++, using ./work/BOUT-dev if available"
  echo "  build_both       Rebuild both hermes and BOUT++"
  echo "  hermes <case>    Run a hermes-3 case on a single rank (e.g. work/test)"
  echo "  hermes_parallel <n> <case>  Run a case on <n> MPI ranks (e.g. 4 work/test)"
  echo "  fix_permissions  Fix ownership of ./work so you can access it"
  echo "  jupyter          Start the Jupyter server on http://localhost:8888"
  echo
  echo "Maintenance:"
  echo "  cleanup          Stop and remove orphaned/stopped containers from this project"
  echo "  rm_docker        Remove ALL containers, volumes and images for this project"
}

# Make sure the environment is set up correctly before running anything.
# Returns non-zero (and explains) if something is missing.
check_environment() {
  if [ ! -f "docker-compose.yaml" ]; then
    warn "Error: no docker-compose.yaml in $PWD."
    warn "Run this from the 'docker' folder of your hermes-3 setup."
    return 1
  fi
  if [ ! -f ".env" ]; then
    warn "Error: no .env file in $PWD."
    warn "Run 'sh setup.sh' (or development-setup.sh) first to generate it."
    return 1
  fi
  if [ ! -d "work" ]; then
    warn "Error: no 'work' subfolder in $PWD to mount into the container."
    warn "Run 'sh setup.sh' (or development-setup.sh) first to create it."
    return 1
  fi
  return 0
}

# Stop and remove orphaned/stopped containers left over from this project.
docker_cleanup() {
  notice "Stopping this project's containers and removing orphans..."
  docker compose down --remove-orphans
  notice "Removing any remaining stopped 'run' containers for this project..."
  docker compose rm --force --stop 2>/dev/null
  notice "Cleanup complete."
}

# Remove everything associated with this project's docker images:
# containers, named/anonymous volumes and the images themselves.
docker_rm() {
  # Image repositories for this project (all tags of these are removed).
  local repos="ghcr.io/boutproject/hermes-3 ghcr.io/boutproject/hermes-3-jupyter"

  warn "This will remove ALL containers, volumes and images (every tag) for this project:"
  for repo in $repos; do
    echo "  - ${repo}"
  done
  printf "Are you sure? [y/N] "
  read -r reply
  case "$reply" in
    [Yy]*) ;;
    *) notice "Aborted."; return 0 ;;
  esac

  notice "Stopping containers and removing volumes/orphans..."
  docker compose down --volumes --remove-orphans --rmi all
  notice "Removing any remaining stopped 'run' containers for this project..."
  docker compose rm --force --stop --volumes 2>/dev/null

  # Remove every tag of these images, in case they weren't created via compose.
  notice "Removing images (all tags)..."
  for repo in $repos; do
    local images
    images=$(docker images --filter "reference=${repo}" --quiet | sort -u)
    if [ -n "$images" ]; then
      # shellcheck disable=SC2086
      docker image rm --force $images 2>/dev/null
    fi
  done

  notice "Purge complete."
  echo
  echo "Note: this removes only this project's images/containers/volumes."
  echo "To reclaim build cache and everything else Docker has accumulated"
  echo "(affects ALL Docker projects on this machine), run:"
  echo "  docker system prune -a --volumes"
}

# Run a hermes-3 case on multiple MPI ranks.
#
# The 'hermes' compose service overrides the image entrypoint, so it does not
# source the spack environment and 'mpirun' is not on its PATH. We therefore run
# mpirun through the 'shell' service (which keeps /entrypoint.sh) instead.
#
# Usage: run_hermes_parallel <n_ranks> <case> [extra BOUT++ args...]
#   e.g. run_hermes_parallel 4 work/fieldline-par-np4
run_hermes_parallel() {
  local nranks="$1"
  local case_dir="$2"

  if [ -z "$nranks" ] || [ -z "$case_dir" ]; then
    warn "Error: hermes_parallel needs a rank count and a case path."
    warn "For example: run_docker hermes_parallel 4 work/test"
    return 1
  fi
  case "$nranks" in
    ''|*[!0-9]*)
      warn "Error: '$nranks' is not a valid number of ranks."
      return 1 ;;
  esac
  # Drop <n_ranks> and <case> so "$@" holds any extra BOUT++ options to forward.
  shift 2

  # Resolve the executable inside the container the same way the image's 'run'
  # command does (prefer the work/ override build if present), then launch it
  # under mpirun. --oversubscribe lets the rank count exceed the host core count.
  # Args are passed positionally ($1=ranks, $2=case, rest=extra) to avoid quoting
  # issues.
  docker compose run --rm shell bash -c '
    nranks="$1"; shift
    case_dir="/hermes_project/$1"; shift
    exe="${HERMES_BUILD_DIR}/hermes-3"
    [ -d "${HERMES_BUILD_DIR_OVERRIDE}" ] && exe="${HERMES_BUILD_DIR_OVERRIDE}/hermes-3"
    if [ ! -x "${exe}" ]; then
      echo "Error: hermes-3 executable not found at ${exe}." >&2
      echo "Build it first, e.g. ./run.sh build_hermes" >&2
      exit 1
    fi
    if [ ! -f "${case_dir}/BOUT.inp" ]; then
      echo "Error: ${case_dir} does not contain a BOUT.inp file." >&2
      exit 1
    fi
    echo "Using ${exe}"
    echo "Running on ${nranks} rank(s) in ${case_dir}"
    exec mpirun --oversubscribe -np "${nranks}" "${exe}" -d "${case_dir}" "$@"
  ' _ "$nranks" "$case_dir" "$@"
}

run_docker() {
  local service="$1"

  if [ -z "$service" ] || [ "$service" = "help" ] || [ "$service" = "-h" ] || [ "$service" = "--help" ]; then
    run_docker_help
    return 0
  fi

  # Maintenance shortcut
  if [ "$service" = "cleanup" ]; then
    check_environment || return 1
    docker_cleanup
    return $?
  fi

  # Maintenance shortcut: remove all images, containers and volumes
  if [ "$service" = "rm_docker" ]; then
    check_environment || return 1
    docker_rm
    return $?
  fi

  # Parallel run: not a compose service, but a helper that drives the 'shell'
  # service so mpirun has the spack environment on its PATH.
  if [ "$service" = "hermes_parallel" ]; then
    check_environment || return 1
    shift
    run_hermes_parallel "$@"
    return $?
  fi

  # Make sure the environment is set up correctly
  check_environment || return 1

  # Make sure the requested action maps to a real compose service
  if ! contains_word "$service" "$VALID_SERVICES"; then
    warn "Error: '$service' is not a valid service."
    echo
    run_docker_help
    return 1
  fi

  # Drop the service name so the rest are arguments passed to the container
  shift

  # Make sure services that need an argument were given one
  if contains_word "$service" "$ARG_REQUIRED_SERVICES" && [ "$#" -eq 0 ]; then
    warn "Error: the '$service' service requires an argument."
    warn "For example: run_docker $service work/test"
    return 1
  fi

  if [ "$service" = "jupyter" ]; then
    # jupyter exposes ports and is meant to stay up
    notice "Starting jupyter (Ctrl-C to stop). Open http://localhost:8888"
    docker compose up jupyter
  else
    # --rm cleans the container up after it exits so nothing is orphaned
    docker compose run --rm "$service" "$@"
  fi
}

# If executed directly (not sourced), run the function with the given arguments.
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  run_docker "$@"
fi
