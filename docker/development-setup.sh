#!/bin/bash
# Sets up a clean development environment for use with the Hermes-3 docker image
# in a "hermes-3-docker" folder. Copy this script to your local machine and execute
# it with "sh development-setup.sh".

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
quiet() {
  "$@" > /dev/null 2>&1
}

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
quiet git remote add -f origin git@github.com:boutproject/hermes-3.git
# We only want to pull the "docker" folder for the project
git config core.sparseCheckout true
echo "docker" >> .git/info/sparse-checkout
quiet git pull origin master
mv docker/* .
rmdir docker
# Uninitialize the git repository
rm -rf .git
sh setup.sh

HERMES_SRC_DIR_OVERRIDE=work/hermes-3
BOUTPP_SRC_DIR_OVERRIDE=work/BOUT-dev
HERMES_CONFIG_OVERRIDE=work/hermes_config.cmake
BOUTPP_CONFIG_OVERRIDE=work/boutpp_config.cmake

notice "Cloning hermes-3/master into $PWD/$HERMES_SRC_DIR_OVERRIDE"
quiet git clone git@github.com:boutproject/hermes-3.git $HERMES_SRC_DIR_OVERRIDE
notice "Copying hermes_config.cmake into $HERMES_CONFIG_OVERRIDE"
cp image_ingredients/hermes_config.cmake $HERMES_CONFIG_OVERRIDE

# Extract the git submodule hash needed by BOUT-dev
BOUT_SUBMODULE_HASH=$(git -C $HERMES_SRC_DIR_OVERRIDE submodule status | grep "BOUT-dev" | awk '{print substr($1, 2)}')

notice "Cloning BOUT-dev/$BOUT_SUBMODULE_HASH into $PWD/$BOUTPP_SRC_DIR_OVERRIDE"
quiet git clone git@github.com:boutproject/BOUT-dev.git $BOUTPP_SRC_DIR_OVERRIDE
quiet git -C $BOUTPP_SRC_DIR_OVERRIDE checkout $BOUT_SUBMODULE_HASH
quiet git -C $BOUTPP_SRC_DIR_OVERRIDE submodule update --init --recursive
notice "Copying boutpp_config.cmake into $BOUTPP_CONFIG_OVERRIDE"
cp image_ingredients/boutpp_config.cmake $BOUTPP_CONFIG_OVERRIDE

if [ -n "$GITHUB_USERNAME" ]; then
    notice "Setting the 'fork' remotes to git@github.com:$GITHUB_USERNAME/hermes-3.git and git@github.com:$GITHUB_USERNAME/BOUT-dev.git."
    git -C $HERMES_SRC_DIR_OVERRIDE remote add fork git@github.com:$GITHUB_USERNAME/hermes-3.git
    git -C $BOUTPP_SRC_DIR_OVERRIDE remote add fork git@github.com:$GITHUB_USERNAME/BOUT-dev.git
    quiet git -C $HERMES_SRC_DIR_OVERRIDE fetch fork
    quiet git -C $BOUTPP_SRC_DIR_OVERRIDE fetch fork
else
    warn "Must set the GITHUB_USERNAME variable to set up the 'fork' remotes."
fi

notice "Finished setting up $PWD/hermes-3-docker"
