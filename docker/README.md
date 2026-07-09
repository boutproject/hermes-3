# Docker image for Hermes

This document describes the Docker setup for building and running the Hermes-3 and BOUT++ plasma simulation frameworks. It outlines the purpose of the `Dockerfile`, the `docker-compose.yml` services, and provides instructions on how to interact with the Docker image.

## Overview

This Docker configuration provides a containerized environment for building and running Hermes-3. It's designed to be flexible, allowing you to build the software and run simulations using files from your local machine, and even modify the source code and build configurations within the container by leveraging the `work` folder.

### Multi-architecture support

The published images (`ghcr.io/boutproject/hermes-3`, `hermes-3-builder`, and `hermes-3-jupyter`) are built for both `linux/amd64` and `linux/arm64`. Docker automatically pulls the variant matching your machine, so Apple-silicon Macs get a native `arm64` image and no longer rely on x86 emulation. The Spack build is pinned to a portable microarchitecture baseline per architecture — `aarch64` on arm64 and `x86_64_v3` on amd64 — rather than the build runner's native (SVE2 / AVX-512) target, so the images also run correctly under emulation and on older CPUs. (See the [Spack binary caches](#spack-binary-caches-and-the-microarchitecture-pin) section for the details of this choice.)

## Quick start for developers

For most developers the fastest way in is `development-setup.sh`. It sets up the
docker tooling and wires up editable Hermes-3 and BOUT++ sources that you can rebuild
inside the container, so it is the recommended starting point. Anything it fetches is
cloned over **HTTPS**, so no SSH keys are required. It behaves differently depending
on where you run it:

### In place — inside an existing checkout (recommended)

Run it from the `docker` folder of a checkout you already have. It reuses that
checkout **directly** as the build source — no second copy is cloned — so
`./run.sh build_*` compiles your actual working tree, branch and all:

```bash
cd docker
sh development-setup.sh
```

It detects the surrounding checkout, then:

1. Runs `setup.sh` to create the `.env` file and the `work` folder (see
   [Running without the helper scripts](#running-without-the-helper-scripts)).
2. Writes a `docker-compose.override.yaml` that bind-mounts the checkout (and its
   `external/BOUT-dev` submodule) into the `build_hermes`, `build_boutpp` and
   `build_both` services, pointing their source overrides at it. Compose merges this
   file automatically. Delete the file to revert to the image's built-in sources.
3. Copies the default CMake config files into `work/` so you can tweak build options.

The build stays out of your source tree: build output goes to `work/`
(out-of-source), and the config files disable the in-source
`git submodule update` that CMake would otherwise run at configure time, so **your
checkout — including `external/BOUT-dev`'s checked-out state — is not modified**. The
build services also run as your user (`UID`/`GID` from `.env`), so nothing they write
to `work/` is root-owned. Because the automatic submodule update is off, the
checkout's submodules must already be populated; if any are missing the script tells
you to run `git submodule update --init --recursive` in your checkout first.

To rebuild both your Hermes-3 and your BOUT-dev in one go, use `./run.sh build_both`.
`./run.sh build_hermes` on its own compiles your Hermes-3 tree against the BOUT++
already built into the image, not your `external/BOUT-dev`.

### Standalone — bootstrap from scratch

Download just this one script and run it in an empty directory. It fetches the
`docker` folder from GitHub into a new `hermes-3-docker` folder and then clones the
sources into `work/`:

```bash
curl -O https://raw.githubusercontent.com/boutproject/hermes-3/master/docker/development-setup.sh
sh development-setup.sh
cd hermes-3-docker
```

Here it will, in addition to running `setup.sh` and copying the config files:

* Clone Hermes-3 (`master`) into `work/hermes-3`, and BOUT++ (pinned to the submodule
  commit Hermes-3 expects) into `work/BOUT-dev`, which then override the versions
  baked into the image (see
  [Overriding Default Builds](#overriding-default-builds-using-the-work-folder)).
* If you export `GITHUB_USERNAME` first, add a `fork` remote to each clone pointing at
  your fork:
    ```bash
    GITHUB_USERNAME=your-username sh development-setup.sh
    ```
  Because the clones use HTTPS you can pull without credentials. To *push* to your
  fork, use a [personal access token](https://docs.github.com/en/authentication), or
  switch that one remote to SSH afterwards, e.g.
  `git -C work/hermes-3 remote set-url fork git@github.com:your-username/hermes-3.git`.

Either way, build and run with `run.sh`, described next.

## Building and running with `run.sh`

`run.sh` is the everyday helper for building and running inside the container. It
wraps the raw `docker compose run --rm <service> ...` commands (which are flexible
but easy to get subtly wrong) and validates your input first: that your environment
is set up (a `docker-compose.yaml`, an `.env`, and a `work` folder all exist — i.e.
that `setup.sh` has run), that the requested service is real, and that services which
need an argument (currently `hermes`) were given one.

Use it either directly or by sourcing it:

```bash
./run.sh <service> [arguments]
# or
source run.sh
run_docker <service> [arguments]
```

Common commands:

```bash
./run.sh shell               # interactive shell in the image
./run.sh build_both          # rebuild hermes and BOUT++ (uses ./work if present)
./run.sh hermes work/case    # run a case on a single process
./run.sh hermes work/case 4  # run a case on 4 MPI ranks
./run.sh jupyter             # start the Jupyter server on http://localhost:8888
./run.sh help                # list the available services
```

### Running in parallel

`./run.sh hermes work/case` runs Hermes-3 on a single process. Pass an optional
rank count as a third argument to run it under MPI instead:

```bash
./run.sh hermes work/case 4   # equivalent to: mpirun -np 4 hermes-3 -d work/case
```

The run is launched with `mpirun`. If the rank count fits within the CPUs Docker
exposes to the container it runs normally; if it exceeds them, the run prints a
warning and adds `--oversubscribe` so it still starts (ranks then time-slice the
same CPUs and run slowly — raise Docker's CPU allocation for better performance).
Make sure the case is decomposed for the number of ranks you request
(`NXPE`/`NYPE` in `BOUT.inp`).

### Tidying up

```bash
./run.sh cleanup     # stop and remove orphaned/stopped containers for this project
./run.sh rm_docker   # remove ALL containers, volumes and images for this project
```

`cleanup` runs `docker compose down --remove-orphans` followed by
`docker compose rm`, to clear out containers left behind by interrupted runs.
`rm_docker` goes further and removes this project's images and named volumes too
(after a confirmation prompt).

## Running without the helper scripts

The helper scripts above are thin wrappers over plain `docker compose`. If you only
want to run pre-built simulations — no local source checkout — you can drive it
directly:

1.  **Navigate to the `docker` directory:**
    ```bash
    cd docker
    ```
2.  **Set up the environment:** `setup.sh` writes a `.env` file (your user/group IDs,
    the image tags and build settings) and creates the `work` folder that is mounted
    into the container. `work` is your central point for moving simulation files to
    and from the container, and for overriding the default builds.
    ```bash
    ./setup.sh
    ```
3.  **Prepare your simulation case:** create a subdirectory of `work` and copy your
    `BOUT.inp` into it:
    ```bash
    mkdir -p work/case
    cp /path/to/your/BOUT.inp work/case/
    ```
4.  **Run the simulation:** this starts a temporary container from the `hermes`
    service, runs Hermes-3 against `work/case`, and removes the container when it
    finishes. Output files are written back to `work/case` on your machine.
    ```bash
    docker compose run --rm hermes work/case      # single process
    docker compose run --rm hermes work/case 4    # 4 MPI ranks
    ```
    These are exactly the commands `./run.sh hermes work/case [n]` runs for you.

## Selecting an image tag

By default the services use the `latest` tag of `ghcr.io/boutproject/hermes-3` and `ghcr.io/boutproject/hermes-3-jupyter`. The tags are controlled by two variables, `HERMES_TAG` (for the `hermes-3` image, used by every service except `jupyter`) and `JUPYTER_TAG` (for the `hermes-3-jupyter` image). `setup.sh` writes both to `.env` set to `latest`.

To use a different tag, either:

* **Edit `.env`** (or re-run `./setup.sh` and then edit it) to set a persistent value:
    ```
    HERMES_TAG=edge
    JUPYTER_TAG=latest
    ```
* **Override on the command line** for a single invocation:
    ```bash
    HERMES_TAG=edge docker compose run --rm hermes work/case
    ```
    Note that a value set in `.env` takes precedence over an exported shell variable of the same name, so pass the override inline (as above) rather than `export`-ing it if you have a value in `.env`.

If neither the variable nor `.env` provides a value, the compose file falls back to `latest`.

## Controlling build parallelism

When you rebuild BOUT++ or Hermes-3 in the container (`build_boutpp`, `build_hermes`,
`build_both`), the number of parallel compile jobs is controlled by the
`HERMES_BUILD_JOBS` variable, which defaults to `4`. `setup.sh` writes
`HERMES_BUILD_JOBS=4` to `.env`. A low default keeps memory use modest —
BOUT++/Hermes-3 are template-heavy, so each job can consume a lot of RAM — but
if you have spare cores and memory you can raise it:

* **Edit `.env`** for a persistent value: `HERMES_BUILD_JOBS=8`
* **Override for a single build**: `HERMES_BUILD_JOBS=8 ./run.sh build_both`

The same variable is also a build-arg for the image build itself
(`docker build --build-arg HERMES_BUILD_JOBS=8 ...`), where it likewise defaults
to `4`.

## Overriding Default Builds using the `work` Folder

The `work` folder on your local machine is mounted to the `/hermes_project/work` directory inside the Docker container. This allows you to override the default builds of BOUT++ and Hermes-3 that are included in the image by placing your own source code and configuration files within the `work` folder. The build scripts inside the container are designed to prioritize these files if they exist.

Here's how you can leverage the `work` folder to customize the build:

1.  **Modifying BOUT++ Source Code:**
    * Clone the BOUT++ repository into the `work` folder on your local machine. **Important:** The directory must be named exactly `BOUT-dev`.
        ```bash
        git clone https://github.com/boutproject/BOUT-dev.git work/BOUT-dev
        ```
    * When you run `docker compose run --rm build_boutpp` or `docker compose run --rm build_both`, the build process will detect the `BOUT-dev` directory in `work` and use this source code instead of the version built into the image.

2.  **Modifying Hermes-3 Source Code:**
    * Clone the Hermes-3 repository into the `work` folder on your local machine. **Important:** The directory must be named exactly `hermes-3`.
        ```bash
        git clone https://github.com/boutproject/hermes-3.git work/hermes-3
        ```
    * When you run `docker compose run --rm build_hermes` or `docker compose run --rm build_both`, the build process will detect the `hermes-3` directory in `work` and use this source code instead of the default.

3.  **Providing Custom Configuration Files:**
    * If you need to change the CMake build options for BOUT++, create a file named `boutpp_config.cmake` and place it directly inside the `work` folder on your local machine.
    * Similarly, for custom Hermes-3 build options, create a file named `hermes_config.cmake` and place it directly inside the `work` folder.
    * When the build scripts run, they will check for these configuration files in the `work` folder and use them in preference to the default configuration files included in the image.

By placing your source code and configuration files in the correctly named locations within the `work` folder, you can effectively override the default build process and use your own customized versions of BOUT++ and Hermes-3 within the Docker environment without needing to rebuild the entire Docker image.

## Interacting with the Image using `image` Commands

The Docker image includes a set of pre-defined commands accessible via the `/bin/image` script within the container. These commands simplify common tasks such as building the code and running simulations. When you use `docker compose run --rm <service> <command> <arguments>`, you are essentially executing this `/bin/image` script with the specified `<command>` and `<arguments>`. The `image` command is defined by the `image_ingredients/docker_image_commands.sh` file.

Here's a breakdown of the available commands:

* **`build_boutpp`**: This command is used to configure and build the BOUT++ framework. It looks for the BOUT++ source code in `/hermes_project/work/BOUT-dev` (if this directory exists; otherwise, it uses the version included when the image was built). It uses the configuration file located at `/hermes_project/work/boutpp_config.cmake` (if it exists; otherwise, it uses the default). The build output is placed in `/hermes_project/build/boutpp-build`.

* **`build_hermes`**: This command configures and builds the Hermes-3 model. It searches for the Hermes-3 source code in `/hermes_project/work/hermes-3` (overriding the built-in version if present) and uses the configuration file at `/hermes_project/work/hermes_config.cmake` (or the default). It depends on a previously built BOUT++ (using the `build_boutpp` command or the BOUT++ built into the image) and places the resulting executable in `/hermes_project/build/hermes-3-build`.

* **`build_both`**: This is a convenience command that sequentially executes `build_boutpp` followed by `build_hermes`.

* **`run <directory> [nprocs]`**: This command executes the compiled Hermes-3 model. You must provide a `<directory>` argument, which should correspond to a subdirectory within the `work` folder on your local machine (e.g., `work/case`). The command checks for a `BOUT.inp` file within this directory and then runs the Hermes-3 executable, passing the specified directory as the working directory for the simulation. An optional second argument sets the number of MPI ranks: with no value (or `1`) it runs on a single process, and with `n > 1` it launches `mpirun -np n` (adding `--oversubscribe`, with a warning, only if `n` exceeds the cores available to the container).

* **`build_hermes_and_run <directory>`**: This command first builds the Hermes-3 model (using the same logic as `build_hermes`) and then immediately runs it in the specified `<directory>` (similar to the `run` command).

* **`build_both_and_run <directory>`**: This command first builds both BOUT++ and Hermes-3 (like the `build_both` command) and then executes the Hermes-3 model in the provided `<directory>`.

* **`fix_permissions`**: This command adjusts the ownership of the `/hermes_project/work` directory within the container to match the user and group IDs of your host machine. This is useful when you encounter permission issues when trying to access or modify files in the `work` folder from within a container started with root privileges (e.g., using `docker compose run --rm sudo`).

When you use commands like `docker compose run --rm hermes work/case`, the `docker-compose.yml` defines the `hermes` service in a way that effectively runs `/bin/image run work/case` inside the container. Similarly, `docker compose run --rm build_hermes` executes `/bin/image build_hermes`.

## Continuous Integration and Delivery

All published images are built and pushed to the GitHub Container Registry (`ghcr.io/boutproject/...`) by GitHub Actions. There are three images and four workflows.

### The three images and how they relate

* **`hermes-3-builder`** (`docker/hermes-3-builder.dockerfile`) — a heavyweight image containing the full Spack-built scientific toolchain (MPI, PETSc, SLEPc, SUNDIALS, netCDF, cmake, Python, …). This is the slow, expensive image to build, so it is built infrequently and cached aggressively (see [Spack binary caches](#spack-binary-caches-and-the-microarchitecture-pin)).
* **`hermes-3`** (`docker/hermes-3.dockerfile`) — the runtime image. It starts `FROM ghcr.io/boutproject/hermes-3-builder:${BUILDER_TAG}` to pull in the toolchain, then compiles BOUT++ and Hermes-3 on top of it and copies the result into a slim `ubuntu:24.04` runtime stage. Because it builds on the builder, the builder image must exist first.
* **`hermes-3-jupyter`** (`docker/hermes-3-jupyter.dockerfile`) — an independent image built `FROM quay.io/jupyter/scipy-notebook` with `xhermes` and docs tooling added. It does not depend on the other two.

### Per-architecture builds and manifest merging

The multi-arch images are **not** built with QEMU emulation for the heavy compilation (that would be prohibitively slow). Instead, `build_builder_image.yml` and `build_docker_image.yml` use a matrix that builds each architecture on a *native* runner — `ubuntu-latest` for `amd64` and `ubuntu-24.04-arm` for `arm64`. Each arch build pushes its image *by digest only* (no tag), and a following `merge` job assembles the per-arch digests into a single multi-arch manifest and applies the human-readable tags (`latest`, `edge`, semver, etc.). Docker then serves the matching architecture automatically on `docker pull`.

The `hermes-3-jupyter` image is the exception: its steps are only `apt`/`pip` installs, so emulation is cheap enough that `build_jupyter_image.yml` just builds `linux/amd64,linux/arm64` together via QEMU in a single job.

All workflows also use a GitHub Actions layer cache (`cache-from`/`cache-to: type=gha`, scoped per arch) so unchanged Docker layers are reused between runs.

### The four workflows and their triggers

| Workflow | Builds | Triggers |
| --- | --- | --- |
| `build_builder_image.yml` | `hermes-3-builder` (→ `latest`/`edge`/date tags), then chains into building `hermes-3` on top | monthly schedule (1st of the month), `release`, `workflow_dispatch` |
| `build_docker_image.yml` | `hermes-3` only, `FROM` the `latest` builder (→ `latest`/`edge`/semver tags) | push to `master`, `release`, `workflow_dispatch` |
| `build_jupyter_image.yml` | `hermes-3-jupyter` | push to `master`, `release`, `workflow_dispatch` |
| `test_docker_build.yml` | builder and/or `hermes-3` under the `experimental_docker_build` tag (validation only) | pull requests touching `docker/**` or that workflow, `workflow_dispatch` |

The split reflects build cost: the expensive builder is only refreshed on a schedule/release, while the cheaper `hermes-3` runtime image is rebuilt on every push to `master` (picking up source changes) `FROM` the most recent `latest` builder.

### The `experimental_docker_build` tag

`test_docker_build.yml` runs on pull requests that touch the `docker/` folder. It validates that the images still build and publishes the results under a dedicated `experimental_docker_build` tag, so the builder→final chain can be exercised end to end without disturbing the `latest`/`edge` tags that real users rely on. A `changes` job (using `dorny/paths-filter`) decides what to rebuild:

* If a change touches the builder image or anything it depends on (the builder Dockerfile, `spack.yaml`, `spack_config.yaml`, the entrypoint, or the workflow), the builder is rebuilt and published as `experimental_docker_build`, and the `hermes-3` image is then built on top of that freshly built builder.
* If only the final image is affected, `hermes-3` is built by itself, `FROM` the `experimental_docker_build` builder published by a previous run.

### Bootstrapping the `experimental_docker_build` builder

The final-image-only path of `test_docker_build.yml` builds `FROM ghcr.io/boutproject/hermes-3-builder:experimental_docker_build`, so that tag must already exist in the registry. The first time you use this workflow (or any time the experimental builder tag has been deleted), you must seed it by running a build that produces the builder image. The reliable way to do this is to open (or push to) a branch whose changes touch a builder dependency — e.g. `docker/hermes-3-builder.dockerfile` or `docker/image_ingredients/spack.yaml` — which makes the workflow rebuild the builder, publish the `experimental_docker_build` builder tag, and then build the final image on top of it.

The workflow can also be started manually via **workflow_dispatch**, but the builder is only rebuilt if a builder dependency differs from the base branch (the change-detection step compares against the base), so a manual run is not by itself guaranteed to seed the tag.

Note that publishing to the registry requires write access, so this workflow only works for branches pushed to this repository — pull requests opened from forks cannot publish the experimental tag.

### Spack binary caches (and the microarchitecture pin)

To avoid recompiling the whole scientific stack on every build, the builder image pulls prebuilt binaries from two mirrors:

* **Spack's public mirror** (`binaries.spack.io/develop`) — covers most common dependencies.
* **A self-hosted OCI cache** (`ghcr.io/boutproject/hermes-3-spack-cache`) — fills in the config-specific packages the public mirror doesn't carry (e.g. `meson`, `py-meson-python`, `py-numpy`). Each builder run pulls from it before building and pushes newly-built specs back, so it warms up over time. This is only used when registry credentials are available (same-repo branches); fork PRs and local builds fall back to the public mirror plus compiling from source.

Both caches are keyed on each package's Spack concretization hash, which **includes the microarchitecture target** — and this is exactly why the target is pinned.

#### Why a hard target pin (not just `granularity: generic`)

The images need to run on older CPUs and under emulation, so they must not contain instructions (SVE2 on arm64, AVX-512 on amd64) that the runtime host lacks. `spack.yaml` sets `concretizer.targets.granularity: generic` as a first line of defence, but **that setting alone is not enough**: it is only a *preference*, and with cache reuse enabled (which the OCI cache relies on) the concretizer will happily reuse a prebuilt binary at a newer native microarchitecture. In practice the public mirror serves arm64 binaries built for `armv9.0a` (SVE2), so a build that merely *preferred* generic still pulled them — and the resulting image crashed with `Illegal instruction` the moment the solver ran on Apple-silicon / a Docker VM.

The fix is a **hard target constraint**, which cache reuse cannot override. Because the shared `spack.yaml` can't branch on architecture, `hermes-3-builder.dockerfile` injects it per-arch before installing:

```sh
case "$(uname -m)" in
  aarch64) HERMES_TARGET=aarch64 ;;
  x86_64)  HERMES_TARGET=x86_64_v3 ;;
esac
spack -e . config add "packages:all:require:target=${HERMES_TARGET}"
```

`aarch64` and `x86_64_v3` are portable baselines (`x86_64_v3` adds AVX2, BMI1/2 and FMA on top of x86-64_v2 ISA). Requiring a *generic-family* target this way does **not** break concretization of the externally-found gcc — a newer compiler can always target an older ISA — despite the compiler being a graph node in Spack 1.x.

Because the target is now fixed regardless of which runner a job lands on, concretization — and therefore OCI-cache hits — are **deterministic across the heterogeneous amd64 runner fleet**. The one residual cost: the public mirror ships some amd64 binaries only at other `x86_64_vN` levels, so a few packages may compile from source at the `v3` floor; they populate the OCI cache on the first run and are reused incrementally thereafter (arm64, already generic on the public mirror, is fully cache-served from the start).

#### Saving build progress on failure

In the OCI-cache branch the builder runs `spack install` **without** `--fail-fast` and pushes to the OCI cache **regardless of the install exit code**, then re-propagates that code so a failed build still fails CI:

```sh
{ spack install ; rc=$? ; \
  spack buildcache push --unsigned --update-index hermes-oci || true ; \
  exit $rc ; }
```

Previously `spack install --fail-fast && spack buildcache push` discarded *every* package built during a run that failed near the end. Now whatever installed is cached, so a retry reuses it as prebuilt binaries and only recompiles what actually failed.

## How the Docker Images are Built (Dockerfiles Explained)

The runtime image is built from **two** Dockerfiles: `docker/hermes-3-builder.dockerfile` produces the `hermes-3-builder` image (the scientific toolchain), and `docker/hermes-3.dockerfile` produces the final `hermes-3` image `FROM` that builder. Splitting them lets the expensive toolchain be built and cached independently of the frequently-changing application build.

### The builder image (`hermes-3-builder.dockerfile`)

1.  **Base image (`FROM spack/ubuntu-noble@sha256:...`)**: a Spack image pinned by digest, based on Ubuntu 24.04. It already ships Spack *and* a full GCC toolchain (`gcc`/`g++`/`gfortran`/`make`, `libc6-dev`) plus `git`, so **no `apt-get` step is needed** — nothing is installed via the OS package manager.
2.  **Spack config + external compiler**: the global Spack config (`spack_config.yaml`, which sets the install tree to `/opt/software`) is copied in, and `spack external find gcc` registers the base image's GCC as an external package. In Spack 1.x compilers are graph nodes and the concretizer only accepts them when found as external *packages*, not via the legacy `spack compiler find`.
3.  **Binary mirrors**: Spack's signed public mirror (`binaries.spack.io/develop`) is added and its keys trusted, and — when registry credentials are supplied — the self-hosted OCI cache is added too. See the [Spack binary caches](#spack-binary-caches-and-the-microarchitecture-pin) section.
4.  **Pin the microarchitecture target**: a hard `packages.all.require: target=…` (per-arch: `aarch64` / `x86_64_v3`) is added to the environment so cache reuse can't pull a non-portable native build. See the [Spack binary caches](#spack-binary-caches-and-the-microarchitecture-pin) section.
5.  **Install the environment (`spack.yaml` → `spack install`)**: the environment manifest (`docker/image_ingredients/spack.yaml`) lists the desired packages (`cmake`, `python`, MPI, PETSc, SLEPc, SUNDIALS, netCDF, …) and the base concretizer settings (`require: %gcc`, `granularity: generic`). Note that **`cmake` is built by Spack** here — it is not in the base image. Whatever installs is pushed back to the OCI cache (even if the install later fails); see [saving build progress on failure](#saving-build-progress-on-failure).
6.  **Activation script (`spack env activate --sh -d . > activate.sh`)** and **entrypoint**: `activate.sh` activates the Spack environment when sourced, and `docker_entrypoint.sh` (which sources it, then `exec "$@"`) is installed as the image entrypoint.

### The final runtime image (`hermes-3.dockerfile`)

1.  **`FROM ghcr.io/boutproject/hermes-3-builder:${BUILDER_TAG}` (as `builder`)**: pulls in the toolchain. `BUILDER_TAG` is a build-arg (default `latest`; CI overrides it, e.g. to `edge` or `experimental_docker_build`).
2.  **`FROM ubuntu:24.04`** — a fresh, slim runtime stage. The Spack environment, `/opt/software`, and `/opt/views` are copied from the `builder` stage, so the built toolchain is available without the Spack machinery.
3.  **Runtime tools (`apt-get ... git vim ca-certificates build-essential`)**: a small set of tools useful for interacting with and rebuilding inside the container. (`cmake` is *not* installed here — it comes from the Spack environment.)
4.  **Working directory and environment variables (`WORKDIR /hermes_project`, `ENV ...`)**: define the locations of the source, build, and config directories for BOUT++ and Hermes-3. The `*_OVERRIDE` variables point into `work/`, so mounted local source and config files take precedence over the built-in versions (see [Overriding Default Builds](#overriding-default-builds-using-the-work-folder)).
5.  **Copy source and init submodules (`COPY . ${HERMES_SRC_DIR}`, `git submodule update ...`)**: the Hermes-3 source is copied in (files excluded via `docker/hermes-3.dockerfile.dockerignore`) and its submodules (including BOUT++) are initialized.
6.  **Copy default config files** (`boutpp_config.cmake`, `hermes_config.cmake`), used unless overridden from `work/`.
7.  **Configure and build BOUT++ then Hermes-3 (`. activate.sh && cmake ... && cmake --build ... --parallel ${HERMES_BUILD_JOBS}`)**: both are compiled at image-build time. `-DCMAKE_PREFIX_PATH` points Hermes-3 at the freshly-built BOUT++. `HERMES_BUILD_JOBS` (a build-arg, default `4`) sets the number of parallel compile jobs.
8.  **Helper commands and entrypoint**: `docker_image_commands.sh` is installed as `/bin/image` (providing `build_boutpp`, `build_hermes`, `run`, etc.), `docker_entrypoint.sh` becomes the entrypoint, and the default command is `/bin/bash`.

## Required upkeep

The Docker images pin their inputs deliberately (for reproducible, cache-friendly
builds), which means they don't update themselves. The following need periodic
manual attention to keep the images building, secure, and current:

* **Spack base image** (`docker/hermes-3-builder.dockerfile`): pinned as
  `spack/ubuntu-noble:<tag>@sha256:<digest>`. Docker resolves purely by digest,
  so bump the tag **and** the digest together by hand. Watch Spack release notes
  when bumping — e.g. 1.1 had a PETSc packaging bug that broke this build and 1.2
  fixed it, so a version jump can require touching `spack.yaml` too.
* **Jupyter base image** (`docker/hermes-3-jupyter.dockerfile`): pinned as
  `quay.io/jupyter/scipy-notebook:<date-tag>@sha256:<digest>` (the Docker Hub
  `jupyter/*` images are frozen at Oct 2023 — stay on quay.io). Bump the dated
  tag and digest together.
* **Runtime base image** (`docker/hermes-3.dockerfile`): `FROM ubuntu:24.04`.
  Move to the next Ubuntu LTS when appropriate, keeping it aligned with the
  builder's Ubuntu release.
* **CI runner architecture labels** (`.github/workflows/build_*image.yml`,
  `test_docker_build.yml`): the per-arch build matrices pin `ubuntu-24.04` and
  `ubuntu-24.04-arm`. GitHub retires old runner images, so bump these to the
  next version when it lands — and keep the amd64 and arm64 labels on the same
  Ubuntu release so the two architectures stay aligned.
* **Microarchitecture target floor** (`hermes-3-builder.dockerfile`, per-arch
  `x86_64_v3` / `aarch64`): revisit if the portability baseline needs to change
  (see [the microarchitecture pin](#why-a-hard-target-pin-not-just-granularity-generic)).
* **Spack environment** (`docker/image_ingredients/spack.yaml`): the package set
  and versions used to build the scientific stack. Refresh as dependencies
  (PETSc, SUNDIALS, MPI, …) move forward.

After changing any pin, the `test_docker_build.yml` workflow will run to confirm everything builds correctly before it reaches `master` and updates the default `latest` image.

## Help!!! I don't have permission to delete the `hermes-3-docker/work` folder

One somewhat annoying problem with using this docker image is the possibility that you end up with files in the `hermes-3-docker/work` folder that you don't have permission to delete. If you still have the `docker-compose.yaml` and `.env` file available, you can run `docker compose run --rm fix_permissions` and then proceed with deleting the `hermes-3-docker/work` folder. However, if you've already deleted these files, run the following command
```
docker run --rm -v "${PWD}/hermes-3-docker/work:/hermes_project/work" -e "PUID=$(id -u)" -e "PGID=$(id -g)" ghcr.io/boutproject/hermes-3 image fix_permissions
```
to adjust the permissions of `./hermes-3-docker/work`. You should then be able to `rm -rf hermes-3-docker/work`.
