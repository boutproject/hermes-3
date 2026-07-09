#!/bin/bash
# Helper for running the Hermes-3 docker compose services.
#
# Usage:
#   Source this file and call run_docker, or run it directly:
#     ./run.sh <service> [arguments]
#
# Examples:
#   ./run.sh shell                 # interactive shell in the image
#   ./run.sh build_both            # rebuild hermes and BOUT++
#   ./run.sh hermes work/test      # run a hermes-3 case on a single process
#   ./run.sh hermes work/test 4    # run a hermes-3 case on 4 MPI ranks
#   ./run.sh jupyter               # start the Jupyter server
#   ./run.sh cleanup               # tidy up orphaned/stopped containers
#   ./run.sh rm_docker             # remove all images, containers and volumes

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
# Services that compile in parallel and honour HERMES_BUILD_JOBS
BUILD_SERVICES="build_hermes build_boutpp build_both"

# Return 0 if $1 appears as a whitespace-separated word in $2
contains_word() {
  case " $2 " in
    *" $1 "*) return 0 ;;
    *) return 1 ;;
  esac
}

# Resolve a variable the way docker compose does: a value exported in the
# environment wins, otherwise fall back to the value assigned in .env. Uses
# printenv (not $VAR) deliberately -- bash always sets UID/EUID as shell
# variables but does not export them, so they are NOT what compose (a child
# process) actually sees; printenv reports only the exported environment.
resolve_env() {
  local key="$1" value
  value=$(printenv "$key" 2>/dev/null)
  if [ -n "$value" ]; then
    printf '%s' "$value"
    return 0
  fi
  sed -n "s/^${key}=//p" .env 2>/dev/null | tail -n1
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
  echo "  hermes <case> [n]  Run a hermes-3 case (path required, e.g. work/test)."
  echo "                     Optional n runs it on n MPI ranks (default 1)."
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
  # Every non-root service runs as ${UID}:${GID} (from .env) so build/run output
  # on the mounted work/ is owned by you. If .env is empty or truncated, compose
  # substitutes a blank user and the container runs as root, littering root-owned
  # files across work/ on the host. Refuse to start unless both are integers.
  local uid gid bad=""
  uid=$(resolve_env UID)
  gid=$(resolve_env GID)
  case "$uid" in "" | *[!0-9]*) bad="UID='${uid}'" ;; esac
  case "$gid" in "" | *[!0-9]*) bad="${bad:+$bad, }GID='${gid}'" ;; esac
  if [ -n "$bad" ]; then
    warn "Error: ${bad} is not a non-negative integer (UID/GID come from .env)."
    warn "Without a numeric UID:GID the container runs as root and creates root-owned files in work/."
    warn "Run 'sh setup.sh' to regenerate .env."
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

# Fail early if the requested parallel build job count exceeds the CPUs Docker
# exposes. This mirrors the container's --oversubscribe check for MPI runs, but
# from the host side. HERMES_BUILD_JOBS is resolved the same way docker compose
# resolves it: the shell environment wins, otherwise the value in .env, otherwise
# 4. The core count comes from `docker info` (the CPUs available to the Docker
# engine/VM, which is what the build container sees). Over-requesting would
# oversubscribe and thrash the CPU, so we refuse it. Returns non-zero to abort.
check_build_jobs() {
  local jobs
  jobs=$(resolve_env HERMES_BUILD_JOBS)
  jobs="${jobs:-4}"

  # A non-integer (or zero) job count is bad config that would fail the build
  # later with a murkier error (cmake --parallel abc), so reject it up front.
  case "$jobs" in
    "" | *[!0-9]* | 0*)
      warn "Error: HERMES_BUILD_JOBS='${jobs}' is not a positive integer."
      warn "Set it to the number of parallel build jobs (e.g. 4), in .env or the environment."
      return 1 ;;
  esac
  # The core count comes from `docker info`; if it isn't available (e.g. the
  # Docker daemon isn't running) we can't compare, so skip rather than block --
  # docker compose will report the real problem.
  local ncores
  ncores=$(docker info --format '{{.NCPU}}' 2>/dev/null)
  case "$ncores" in
    "" | *[!0-9]*) return 0 ;;
  esac

  if [ "$jobs" -gt "$ncores" ]; then
    warn "Error: HERMES_BUILD_JOBS=${jobs} exceeds the ${ncores} CPU(s) Docker exposes."
    warn "Oversubscribed build jobs thrash the CPU. Lower HERMES_BUILD_JOBS to at most ${ncores}, or raise Docker's CPU allocation."
    return 1
  fi
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

  # For 'hermes', an optional second argument is the MPI rank count; if given it
  # must be a positive integer (matches the container's ^[1-9][0-9]*$ check, so
  # bad input fails here before a container is even started).
  if [ "$service" = "hermes" ] && [ -n "$2" ]; then
    case "$2" in
      *[!0-9]* | "" | 0*) warn "Error: rank count must be a positive integer, got '$2'."
                          warn "For example: run_docker hermes work/test 4"
                          return 1 ;;
    esac
  fi

  # For build services, abort if the parallel job count exceeds Docker's CPUs.
  if contains_word "$service" "$BUILD_SERVICES"; then
    check_build_jobs || return 1
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
