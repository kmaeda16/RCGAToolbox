# libSBML 5.18.0 MATLAB binaries for Apple Silicon

This directory contains files used to build native Apple Silicon
(`maca64`) MATLAB MEX binaries for libSBML 5.18.0.

## Included files

- `libsbml-5.18.0-maca64.patch`
  Compatibility patch for Apple Silicon and modern CMake.

- `build-libsbml-maca64.sh`
  Build script for generating:
  - `OutputSBML.mexmaca64`
  - `TranslateSBML.mexmaca64`

- `OutputSBML.mexmaca64`
- `TranslateSBML.mexmaca64`

## Source

The libSBML source code itself is not included in this directory.

Download and extract the official libSBML 5.18.0 source distribution from the libSBML project website:

https://sourceforge.net/projects/sbml/files/libsbml/5.18.0/stable/

Use the source archive:

libSBML-5.18.0-core-plus-packages-src.tar.gz

After extraction, the source directory should look like:

libSBML-5.18.0-Source/

## Build requirements

- Apple Silicon Mac
- Apple Silicon version of MATLAB
- Xcode / Apple Clang
- CMake

## Build

Run:

./build-libsbml-maca64.sh \
    /path/to/libSBML-5.18.0-Source \
    /Applications/MATLAB_R20XXx.app

The script:

1. applies the Apple Silicon compatibility patch;
2. configures libSBML with MATLAB support;
3. enables the FBC, Groups, and Qual packages;
4. builds native arm64 MATLAB MEX binaries;
5. copies the resulting binaries into this directory.

## Package configuration

The binaries are built with:

- FBC
- Groups
- Qual

This matches the package configuration of the official
libSBML 5.18.0 MATLAB binaries for Intel macOS.

The configuration can be checked in MATLAB using:

info = OutputSBML()

Expected output includes:

    libSBML_version_string: '5.18.0'
    isFBCEnabled: 'enabled'
    packagesEnabled: 'fbc;groups;qual'

The XML parser version may differ depending on the build environment.

## Tested environment

Tested with:

- Apple Silicon Mac
- MATLAB R2026b
- Xcode 26.6
- libSBML 5.18.0

## Notes

These binaries are intended for use by RCGAToolbox on Apple Silicon.
Most RCGAToolbox users should not need to rebuild them manually.
