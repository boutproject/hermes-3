# Docker image for Hermes

This document describes the Docker setup for building and running the Hermes-3 and BOUT++ plasma simulation frameworks. It outlines the purpose of the `Dockerfile`, the `docker-compose.yml` services, and provides instructions on how to interact with the Docker image.

## Overview

This Docker configuration provides a containerized environment for building and running Hermes-3. It's designed to be flexible, allowing you to build the software and run simulations using files from your local machine.

## Quickstart

### Running a simulation

The fastest way to run a Hermes-3 simulation is
1. Run `./setup.sh` to make a `.env` file with environment variables needed by `docker` and make a `work` subfolder in your current working directory.
2. Copy a `BOUT.inp` file into `work/case` to configure a simulation.
3. Run `docker compose run --rm hermes work/case` to run a Hermes-3 case in the `work` folder

### Modifying the source code of Hermes-3

If you want to directly edit the source code of `hermes-3`
1. Run `./setup.sh` to make a `.env` file with environment variables needed by `docker` and make a `work` subfolder in your current working directory.
2. Get a copy of the source code using `https://github.com/boutproject/hermes-3.git work/hermes-3`. Note: you must call your source code `work/hermes-3`, otherwise it will be ignored.
3. Run `docker compose run --rm build_hermes` to rebuild `hermes-3`. Note: this will use the version of `BOUT-dev` which was used by `hermes-3` when the docker image was built.
4. Follow the steps from the 'running a simulation' section to use the updated executable.

### Modifying the source code of BOUT++

If you want to directly edit the source code of `BOUT-dev`
1. Run `./setup.sh` to make a `.env` file with environment variables needed by `docker` and make a `work` subfolder in your current working directory.
2. Get a copy of the source code using `https://github.com/boutproject/BOUT-dev.git work/BOUT-dev`. Note: you must call your source code `work/BOUT-dev`, otherwise it will be ignored.
3. Run `docker compose run --rm build_both` to rebuild `BOUT-dev` and `hermes-3`.
4. Follow the steps from the 'running a simulation' section to use the updated executable.

### Interact with the image via a terminal

If you'd like fine-tuned control of the image, you can start an interactive terminal using either
* `docker compose run --rm shell` to start a session with user and group id matching your host machine
* `docker compose run --rm sudo` to start a session with root access. Note that this might lead to permissions errors when trying to interact with the `work` folder on your host machine. You can fix these using `image fix_permissions` from a `sudo` terminal session, or `docker compose run --rm fix_permissions` from your host machine.

### Tidying up

In case you want to tidy up after using the image, you can use `docker compose down --remove-orphans` to remove instances which are no longer being used.

You can also
* `docker system prune --volumes` to remove all unused containers, networks and images.

### Building an image locally

To build the `hermes3` image locally, you can use the command `sudo docker build -t hermes3 docker` from the `hermes-3` directory. Note: the `spack install` step takes several hours to run. 
