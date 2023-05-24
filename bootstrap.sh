#!/bin/bash -eu

USAGE="Usage: $0 [options]

Bootstraps the riesling build system (specifies the toolchain for CMake)

Options:
  -f FILE : Use a set of flags for dependency compilation. Currently provided
            options are avx2, abi (for pre-C++11 ABI), and native
"

git submodule update --init --recursive

FLAGS="base"
while getopts "f:hj:" opt; do
    case $opt in
        f) export FLAGS="$OPTARG";;
        h) echo "$USAGE"
           return;;
        j) export VCPKG_MAX_CONCURRENCY="$OPTARG";;
    esac
done
shift $((OPTIND - 1))

# Use Ninja if available, otherwise CMake default
if [ -x "$( command -v ninja )" ]; then
  GEN="-GNinja"
else
  GEN=""
fi

mkdir -p build
cd build
cmake -S ../ $GEN \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE="cmake/toolchain.cmake" \
  -DFLAGS_FILE="${FLAGS}"
cmake --build .
