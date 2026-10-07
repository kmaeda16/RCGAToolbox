#!/bin/bash
set -euo pipefail

# ----------------------------------------------------------------------
# Build libSBML 5.18.0 MATLAB MEX binaries for Apple Silicon (maca64)
#
# Usage:
#   ./build-libsbml-maca64.sh \
#       /path/to/libSBML-5.18.0-Source \
#       /Applications/MATLAB_R2026b.app
#
# The resulting files are copied to the directory containing this script:
#   OutputSBML.mexmaca64
#   TranslateSBML.mexmaca64
# ----------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATCH_FILE="${SCRIPT_DIR}/libsbml-5.18.0-maca64.patch"

if [ "$#" -ne 2 ]; then
    echo "Usage:"
    echo "  $0 /path/to/libSBML-5.18.0-Source /Applications/MATLAB_R20XXx.app"
    exit 1
fi

SOURCE_DIR="$(cd "$1" && pwd)"
MATLAB_ROOT="$2"
BUILD_DIR="${SOURCE_DIR}/build-maca64"

# ----------------------------------------------------------------------
# Check environment
# ----------------------------------------------------------------------

if [ "$(uname -s)" != "Darwin" ]; then
    echo "Error: This script is intended for macOS."
    exit 1
fi

if [ "$(uname -m)" != "arm64" ]; then
    echo "Error: This script must be run natively on Apple Silicon (arm64)."
    exit 1
fi

for cmd in cmake patch file; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Error: '$cmd' was not found."
        exit 1
    fi
done

if [ ! -f "${SOURCE_DIR}/CMakeLists.txt" ]; then
    echo "Error: CMakeLists.txt was not found in:"
    echo "  ${SOURCE_DIR}"
    exit 1
fi

if [ ! -f "${PATCH_FILE}" ]; then
    echo "Error: Patch file was not found:"
    echo "  ${PATCH_FILE}"
    exit 1
fi

if [ ! -f "${MATLAB_ROOT}/bin/maca64/libmex.dylib" ]; then
    echo "Error: Apple Silicon MATLAB was not found at:"
    echo "  ${MATLAB_ROOT}"
    exit 1
fi

echo "libSBML source : ${SOURCE_DIR}"
echo "MATLAB         : ${MATLAB_ROOT}"
echo "Build directory: ${BUILD_DIR}"
echo

# ----------------------------------------------------------------------
# Apply Apple Silicon compatibility patch
# ----------------------------------------------------------------------

cd "${SOURCE_DIR}"

echo "Checking patch ..."

if patch --dry-run -p1 < "${PATCH_FILE}" >/dev/null 2>&1; then
    echo "Applying Apple Silicon patch ..."
    patch -p1 < "${PATCH_FILE}"
elif patch --dry-run -R -p1 < "${PATCH_FILE}" >/dev/null 2>&1; then
    echo "Patch is already applied."
else
    echo "Error: The patch cannot be applied cleanly."
    echo "Make sure that this is an unmodified libSBML 5.18.0 source tree."
    exit 1
fi

# ----------------------------------------------------------------------
# Configure
#
# FBC, Groups, and Qual are enabled to match the official
# libSBML 5.18.0 MATLAB binaries.
# ----------------------------------------------------------------------

echo
echo "Configuring libSBML ..."

rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

cmake .. \
    -DWITH_MATLAB=ON \
    -DENABLE_FBC=ON \
    -DENABLE_GROUPS=ON \
    -DENABLE_QUAL=ON \
    -DMATLAB_ROOT_PATH="${MATLAB_ROOT}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_OSX_ARCHITECTURES=arm64 \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5

# ----------------------------------------------------------------------
# Build
# ----------------------------------------------------------------------

echo
echo "Building libSBML MATLAB bindings ..."

cmake --build . -j 4

OUTPUT_SBML="${BUILD_DIR}/src/bindings/matlab/OutputSBML.mexmaca64"
TRANSLATE_SBML="${BUILD_DIR}/src/bindings/matlab/TranslateSBML.mexmaca64"

if [ ! -f "${OUTPUT_SBML}" ]; then
    echo "Error: OutputSBML.mexmaca64 was not generated."
    exit 1
fi

if [ ! -f "${TRANSLATE_SBML}" ]; then
    echo "Error: TranslateSBML.mexmaca64 was not generated."
    exit 1
fi

# ----------------------------------------------------------------------
# Verify architecture
# ----------------------------------------------------------------------

echo
echo "Checking generated binaries ..."

file "${OUTPUT_SBML}"
file "${TRANSLATE_SBML}"

if ! file "${OUTPUT_SBML}" | grep -q "arm64"; then
    echo "Error: OutputSBML.mexmaca64 is not an arm64 binary."
    exit 1
fi

if ! file "${TRANSLATE_SBML}" | grep -q "arm64"; then
    echo "Error: TranslateSBML.mexmaca64 is not an arm64 binary."
    exit 1
fi

# ----------------------------------------------------------------------
# Copy binaries next to this script
# ----------------------------------------------------------------------

echo
echo "Copying binaries to:"
echo "  ${SCRIPT_DIR}"

cp "${OUTPUT_SBML}" "${SCRIPT_DIR}/"
cp "${TRANSLATE_SBML}" "${SCRIPT_DIR}/"

echo
echo "Build completed successfully."
echo
echo "Generated files:"
echo "  ${SCRIPT_DIR}/OutputSBML.mexmaca64"
echo "  ${SCRIPT_DIR}/TranslateSBML.mexmaca64"
echo
echo "MATLAB verification:"
echo "  info = OutputSBML()"
echo
echo "Expected package configuration:"
echo "  packagesEnabled: 'fbc;groups;qual'"
