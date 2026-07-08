#!/bin/sh
# Seed GNU readline patches into a local Spack source mirror, since Spack fetches
# them from ftpmirror.gnu.org, which is unreachable from our build network.
# Files are laid out by sha256 in Spack's source-cache format; Spack verifies
# the sha256 on use, so transport integrity does not matter.
#
# Usage: seed_gnu_mirror.sh "<GNU_BASE_URL> [GNU_BASE_URL ...]" <MIRROR_DIR> [VERSIONS]
#   GNU_BASE_URLs  space-separated bases serving readline/..., tried in order
#                  per file (e.g. "https://ftp.gnu.org/gnu https://mirrors.kernel.org/gnu")
#   MIRROR_DIR     output dir; register with `spack mirror add <name> file://<dir>`
#   VERSIONS       optional space-separated readline versions to seed patches for
#                  (defaults below); extend when Spack pins a newer readline, or
#                  patches for that version silently won't be seeded.
set -eu

BASES="$1"
MIRROR="$2"
VERSIONS="${3:-8.0 8.1 8.2 8.3 8.4 8.5 8.6}"
ARCHIVE="${MIRROR}/_source-cache/archive"
mkdir -p "${ARCHIVE}"

# Fetch a path (relative to a base) trying each base in turn; on success writes
# to the given output file and returns 0, else returns 1.
fetch_first() {
  _path="$1"; _out="$2"
  for _base in ${BASES}; do
    if curl -fsSL --retry 3 -o "${_out}" "${_base}/${_path}" 2>/dev/null; then
      return 0
    fi
  done
  return 1
}

seeded=0
# readline patches are published as contiguous series per version
# (readline-<v>-patches/readline<vj>-001, -002, ...). Walk each version, and
# within it sequential patch numbers until the first gap (all mirrors miss).
for v in ${VERSIONS}; do
  vj=$(printf '%s' "$v" | tr -d '.')
  n=1
  while :; do
    num=$(printf '%03d' "$n")
    tmp=$(mktemp)
    if fetch_first "readline/readline-${v}-patches/readline${vj}-${num}" "${tmp}"; then
      sha=$(sha256sum "${tmp}" | awk '{print $1}')
      dest="${ARCHIVE}/$(printf '%s' "${sha}" | cut -c1-2)"
      mkdir -p "${dest}"
      mv "${tmp}" "${dest}/${sha}"
      seeded=$((seeded + 1))
      n=$((n + 1))
    else
      rm -f "${tmp}"
      break
    fi
  done
done

echo "seed_gnu_mirror: seeded ${seeded} GNU patch file(s) into ${MIRROR} (mirrors: ${BASES})"
if [ "${seeded}" -eq 0 ]; then
  echo "seed_gnu_mirror: ERROR: nothing seeded; none of the mirrors reachable: ${BASES}" >&2
  exit 1
fi
