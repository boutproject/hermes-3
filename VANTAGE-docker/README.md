# Docker image for hermes-3+vantagereactions

The file `hermes-3-VANTAGE-dependencies.dockerfile` contains the commands to create a docker image which contains the spack dependencies needed to quickly compile hermes-3, BOUT++, NESO-Particles, and VANTAGE-Reactions.

# Make a local docker image

Check for existing docker images.
```
$ docker image ls
```

Make a new docker image from Dockerfile.
```
$ docker build -t hermes_3_vantage_deps_img -f ./path/to/VANTAGE-docker/hermes-3-VANTAGE-dependencies.dockerfile .
```

# Make a local container and test hermes-3 within the container

Check for existing containers
```
$ docker ps -as
```

Make a new container
```
$ cd /path/to/your/hermes-3/
$ docker run --name hermes_3_vantage -v "$(pwd):/root/Hermes-3" -it hermes_3_vantage_deps_img
```

The `-v` option mounts the `hermes-3/` repo to the container on the path `/root/Hermes-3`.

Activate spack.
```
$ . /opt/spack/share/spack/setup-env.sh
```
Activate the environment.
```
$ cd ~/Hermes-3
$ . activate_h3env
```
Install the top-level packages (if spack needs to install more than the top-level, `VANTAGE-Reactions`, `NESO-Particles`, `BOUT++`, `hermes-3` etc, then there is likely a problem with the image, which should capture all build and link dependencies).
```
$ spack concretize -f
$ spack install -j 1
```
Change to the Hermes build directory and test.
```
$ spack cd -b hermes-3
$ ./hermes_unit_tests
$ ctest -j 1
```
To later restart the container.
```
$ docker start -i hermes_3_vantage
```