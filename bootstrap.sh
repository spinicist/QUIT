#!/bin/bash -eux

# This script will configure and build the project using vcpkg for
# dependency management. If you want to build with sanitizer support
# then change the CMAKE_TOOLCHAIN_FILE below to point to the relevant
# .toolchain.cmake file in the cmake/ submodule (or provide your own)

git submodule update --init --recursive
cmake -B build -S . \
    -GNinja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE=cmake/x64-osx-native.toolchain.cmake
cmake --build build
