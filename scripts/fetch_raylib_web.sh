#!/usr/bin/env sh

set -eu

RAYLIB_VERSION="5.5"
RAYLIB_SHA256="aea98ecf5bc5c5e0b789a76de0083a21a70457050ea4cc2aec7566935f5e258e"
RAYLIB_ARCHIVE="raylib-${RAYLIB_VERSION}.tar.gz"
RAYLIB_URL="https://github.com/raysan5/raylib/archive/refs/tags/${RAYLIB_VERSION}.tar.gz"

SCRIPT_DIR=$(CDPATH= cd -P "$(dirname "$0")" && pwd)
PROJECT_ROOT=$(CDPATH= cd -P "${SCRIPT_DIR}/.." && pwd)
DESTINATION=${1:-"${PROJECT_ROOT}/build/web/raylib-${RAYLIB_VERSION}"}

log() {
    printf '%s\n' "$*" >&2
}

fail() {
    log "fetch_raylib_web.sh: $*"
    exit 1
}

is_valid_source_tree() {
    candidate=$1
    [ -f "${candidate}/src/raylib.h" ] &&
        [ -f "${candidate}/src/Makefile" ] &&
        [ -f "${candidate}/LICENSE" ] &&
        grep -Eq '^#define[[:space:]]+RAYLIB_VERSION[[:space:]]+"5\.5"' \
            "${candidate}/src/raylib.h"
}

destination_parent=$(dirname "${DESTINATION}")
destination_name=$(basename "${DESTINATION}")
[ -n "${destination_name}" ] || fail "destination must name a directory"
[ "${destination_name}" != "." ] || fail "destination must name a directory"
[ "${destination_name}" != ".." ] || fail "destination must name a directory"
[ "${destination_name}" != "/" ] || fail "destination must name a directory"

mkdir -p "${destination_parent}"
destination_parent=$(CDPATH= cd -P "${destination_parent}" && pwd)
DESTINATION="${destination_parent}/${destination_name}"

if [ -e "${DESTINATION}" ]; then
    if is_valid_source_tree "${DESTINATION}"; then
        log "Using existing raylib ${RAYLIB_VERSION} source at ${DESTINATION}"
        printf '%s\n' "${DESTINATION}"
        exit 0
    fi
    fail "${DESTINATION} already exists but is not a complete raylib ${RAYLIB_VERSION} source tree"
fi

work_dir=$(mktemp -d "${destination_parent}/.raylib-fetch.XXXXXX") ||
    fail "could not create a temporary directory in ${destination_parent}"

cleanup() {
    if [ -n "${work_dir:-}" ] && [ -d "${work_dir}" ]; then
        rm -rf "${work_dir}"
    fi
}
trap cleanup 0 1 2 15

archive_path="${work_dir}/${RAYLIB_ARCHIVE}"
log "Downloading raylib ${RAYLIB_VERSION} from ${RAYLIB_URL}"

if command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 3 --output "${archive_path}" "${RAYLIB_URL}"
elif command -v wget >/dev/null 2>&1; then
    wget -O "${archive_path}" "${RAYLIB_URL}"
else
    fail "curl or wget is required to download raylib"
fi

if command -v sha256sum >/dev/null 2>&1; then
    actual_sha256=$(sha256sum "${archive_path}" | awk '{print $1}')
elif command -v shasum >/dev/null 2>&1; then
    actual_sha256=$(shasum -a 256 "${archive_path}" | awk '{print $1}')
elif command -v openssl >/dev/null 2>&1; then
    actual_sha256=$(openssl dgst -sha256 "${archive_path}" | awk '{print $NF}')
else
    fail "sha256sum, shasum, or openssl is required to verify raylib"
fi

[ "${actual_sha256}" = "${RAYLIB_SHA256}" ] ||
    fail "SHA-256 mismatch for ${RAYLIB_ARCHIVE}: expected ${RAYLIB_SHA256}, got ${actual_sha256}"

log "Verified ${RAYLIB_ARCHIVE} (${RAYLIB_SHA256})"
tar -xzf "${archive_path}" -C "${work_dir}"

extracted_source="${work_dir}/raylib-${RAYLIB_VERSION}"
is_valid_source_tree "${extracted_source}" ||
    fail "the verified archive did not contain the expected raylib ${RAYLIB_VERSION} source tree"

mv "${extracted_source}" "${DESTINATION}"
log "Installed raylib ${RAYLIB_VERSION} source at ${DESTINATION}"
printf '%s\n' "${DESTINATION}"
