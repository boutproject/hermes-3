#!/bin/bash
# Sets up a clean development environment for use with the Hermes-3 docker image.
#
# Two ways to run it:
#   * Standalone: copy this script to your local machine and run "sh development-setup.sh".
#     It fetches the "docker" folder from GitHub into a new "hermes-3-docker" folder,
#     then clones hermes-3 and BOUT-dev into "work/" as the editable build sources.
#   * In place: run it from an existing checkout's "docker" folder. It skips the fetch
#     and reuses the surrounding checkout (and its external/BOUT-dev submodule) as the
#     build source by generating a docker-compose.override.yaml, so `./run.sh build_*`
#     compiles your working tree directly instead of a separate clone under "work/".

# Define color codes for pretty output
LIGHTRED='\033[1;31m'
LIGHTGREEN='\033[1;32m'
NOCOLOR='\033[0m' # Reset to no color

# Define function for printing messages in different colors
print_message() {
    local color="$1"
    local message="$2"
    printf "${color}%s${NOCOLOR}\n" "$message"
}

# Simplified function calls for different message types
warn() { print_message "${LIGHTRED}" "$1"; }
notice() { print_message "${LIGHTGREEN}" "$1"; }

verbose=false
for arg in "$@"; do
  if [ "$arg" = "-v" ]; then
    verbose=true
  fi
done

# Compact quiet function
quiet() {
  if $verbose; then
    "$@"
  else
    "$@" > /dev/null 2>&1
  fi
}

# Detect whether this script is already being run from inside the "docker"
# folder of a checked-out hermes-3 repo. If so, all the files that the bootstrap
# below would fetch from GitHub are already present, so we skip the fetch and
# just set up the "work" subfolder in place.
in_place=false
if [ -f "docker-compose.yaml" ] && [ -f "setup.sh" ] && [ -d "image_ingredients" ]; then
  in_place=true
  notice "Detected an existing docker folder in $PWD; setting up the 'work' subfolder in place."
else
  if [ -d "hermes-3-docker" ]; then
    warn "Error: Directory '$PWD/hermes-3-docker' already exists. To continue, run"
    warn "cd hermes-3-docker && docker compose run --rm fix_permissions && cd .."
    warn "rm -rf $PWD/hermes-3-docker"
    exit 1
  else
    notice "Setting up a development environment for hermes-3 in $PWD/hermes-3-docker"
  fi
  mkdir -p hermes-3-docker
  cd hermes-3-docker
  quiet git init
  quiet git remote add -f origin https://github.com/boutproject/hermes-3.git
  # We only want to pull the "docker" folder for the project
  git config core.sparseCheckout true
  echo "docker" >> .git/info/sparse-checkout
  quiet git pull origin master
  mv docker/* .
  rmdir docker
  # Uninitialize the git repository
  rm -rf .git
fi
sh setup.sh

# Config files are copied into work/ in both modes so build options can be
# tweaked without touching the source tree.
HERMES_CONFIG_OVERRIDE=work/hermes_config.cmake
BOUTPP_CONFIG_OVERRIDE=work/boutpp_config.cmake

# Decide where the build sources come from. When running in place inside the
# "docker" folder of a real hermes-3 checkout (with the BOUT-dev submodule
# populated), reuse that checkout directly instead of cloning fresh copies into
# work/. Otherwise (standalone bootstrap) clone into work/ as before. The
# in_place guard matters: the standalone path cd's into hermes-3-docker/, whose
# parent might coincidentally be a checkout, and we must not reuse it there.
REPO_ROOT=$(cd .. 2>/dev/null && pwd)
if $in_place && [ -n "$REPO_ROOT" ] && [ -f "$REPO_ROOT/hermes-3.cxx" ] && [ -d "$REPO_ROOT/external/BOUT-dev" ]; then
  notice "Detected the surrounding hermes-3 checkout at $REPO_ROOT."
  notice "Reusing it (and its external/BOUT-dev submodule) as the build source; nothing will be cloned into work/."

  cp image_ingredients/hermes_config.cmake "$HERMES_CONFIG_OVERRIDE"
  cp image_ingredients/boutpp_config.cmake "$BOUTPP_CONFIG_OVERRIDE"

  # These config files set HERMES_UPDATE_GIT_SUBMODULE / BOUT_UPDATE_GIT_SUBMODULE
  # OFF, so the build does NOT run `git submodule update` inside the mounted
  # checkout (it would otherwise reset a work-in-progress external/BOUT-dev). The
  # build is out-of-source (build dir under work/), so the checkout is not written
  # to. But that also means the checkout's submodules must already be populated;
  # warn if any (recursively) are not.
  if git -C "$REPO_ROOT" submodule status --recursive 2>/dev/null | grep -q '^-'; then
    warn "Warning: some submodules in $REPO_ROOT are not initialized."
    warn "Run this in your checkout before building:"
    warn "  git -C $REPO_ROOT submodule update --init --recursive"
  fi

  # Generate a docker-compose override that bind-mounts the checkout into the
  # build services and points the source overrides at it. Compose merges this
  # file automatically, so `./run.sh build_both` compiles your working tree.
  # Build output goes to work/ (out-of-source) and the config files above keep
  # the build from touching the checkout's git state, so the checkout stays
  # clean. Delete this file to fall back to the image's built-in sources.
  # run/shell/jupyter are intentionally left untouched. The volume source is
  # written in long form so paths with spaces or colons are not misparsed by
  # Compose's short "source:target" syntax.
  notice "Writing docker-compose.override.yaml (mounts $REPO_ROOT into the build services)."
  cat > docker-compose.override.yaml <<EOF
# Auto-generated by development-setup.sh for in-place development.
# Bind-mounts the surrounding hermes-3 checkout so the build services compile
# your working tree (and its external/BOUT-dev submodule) instead of copies
# under work/. Delete this file to revert to the image's built-in sources.
services:
  build_hermes:
    volumes:
      - type: bind
        source: "$REPO_ROOT"
        target: /hermes_project/live
    environment:
      - HERMES_SRC_DIR_OVERRIDE=/hermes_project/live
  build_boutpp:
    volumes:
      - type: bind
        source: "$REPO_ROOT"
        target: /hermes_project/live
    environment:
      - BOUTPP_SRC_DIR_OVERRIDE=/hermes_project/live/external/BOUT-dev
  build_both:
    volumes:
      - type: bind
        source: "$REPO_ROOT"
        target: /hermes_project/live
    environment:
      - HERMES_SRC_DIR_OVERRIDE=/hermes_project/live
      - BOUTPP_SRC_DIR_OVERRIDE=/hermes_project/live/external/BOUT-dev
EOF

  notice "Finished setting up $PWD"
  exit 0
fi

# Standalone bootstrap: clone the sources into work/ as separate checkouts.
HERMES_SRC_DIR_OVERRIDE=work/hermes-3
BOUTPP_SRC_DIR_OVERRIDE=work/BOUT-dev

notice "Cloning hermes-3/master into $PWD/$HERMES_SRC_DIR_OVERRIDE"
quiet git clone https://github.com/boutproject/hermes-3.git $HERMES_SRC_DIR_OVERRIDE
notice "Copying hermes_config.cmake into $HERMES_CONFIG_OVERRIDE"
cp image_ingredients/hermes_config.cmake $HERMES_CONFIG_OVERRIDE

# Extract the git submodule hash needed by BOUT-dev
BOUT_SUBMODULE_HASH=$(git -C $HERMES_SRC_DIR_OVERRIDE submodule status | grep "BOUT-dev" | awk '{print substr($1, 2)}')

notice "Cloning BOUT-dev/$BOUT_SUBMODULE_HASH into $PWD/$BOUTPP_SRC_DIR_OVERRIDE"
quiet git clone https://github.com/boutproject/BOUT-dev.git $BOUTPP_SRC_DIR_OVERRIDE
quiet git -C $BOUTPP_SRC_DIR_OVERRIDE checkout $BOUT_SUBMODULE_HASH
quiet git -C $BOUTPP_SRC_DIR_OVERRIDE submodule update --init --recursive
notice "Copying boutpp_config.cmake into $BOUTPP_CONFIG_OVERRIDE"
cp image_ingredients/boutpp_config.cmake $BOUTPP_CONFIG_OVERRIDE

if [ -n "$GITHUB_USERNAME" ]; then
    notice "Setting the 'fork' remotes to https://github.com/$GITHUB_USERNAME/hermes-3.git and https://github.com/$GITHUB_USERNAME/BOUT-dev.git."
    git -C $HERMES_SRC_DIR_OVERRIDE remote add fork https://github.com/$GITHUB_USERNAME/hermes-3.git
    git -C $BOUTPP_SRC_DIR_OVERRIDE remote add fork https://github.com/$GITHUB_USERNAME/BOUT-dev.git
    quiet git -C $HERMES_SRC_DIR_OVERRIDE fetch fork
    quiet git -C $BOUTPP_SRC_DIR_OVERRIDE fetch fork
else
    warn "Must set the GITHUB_USERNAME variable to set up the 'fork' remotes."
fi

notice "Finished setting up $PWD"
