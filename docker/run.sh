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
#   ./run.sh hermes work/test      # run a hermes-3 case (argument required)
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
  echo "  hermes <case>    Run a hermes-3 case (requires a path, e.g. work/test)"
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
