#!/bin/bash -eu

USAGE="Usage: $0 [options]

Bootstraps the QUIT build system (specifies the toolchain for CMake)

Options:
  -i DIR Install riesling to this directory, e.g. $HOME/.local
  -j N   Restrict parallel build to this many threads
  -h     Print this message
"
PAR=""
PREFIX=""
ALL=""
DEBUG="-g"
while getopts "f:hi:j:w" opt; do
    case $opt in
        i) PREFIX="-DCMAKE_INSTALL_PREFIX=$OPTARG -DCMAKE_PREFIX_PATH=$OPTARG";;
        j) export VCPKG_MAX_CONCURRENCY=$OPTARG
           PAR="-j $OPTARG";;
        h) echo "$USAGE"
           return;;
        w) ALL="-Wall"
           DEBUG="-g -fsanitize=address,undefined";;
    esac
done
shift $((OPTIND - 1))

# Use Ninja if available, otherwise CMake default
if [ -x "$( command -v ninja )" ]; then
  GEN="-GNinja"
else
  GEN=""
fi

# If vcpkg is not installed, install it
if [[ (-x "$( command -v vcpkg )") && (-n "${VCPKG_ROOT:-}") ]]; then
  echo "vcpkg installed"
else
  git clone https://github.com/microsoft/vcpkg.git .vcpkg
  cd .vcpkg && ./bootstrap-vcpkg.sh && cd ..
  export VCPKG_ROOT="$PWD/.vcpkg"
  export PATH="$VCPKG_ROOT:$PATH"
fi

mkdir -p build
cmake -S . -B build $GEN \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE="${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_INSTALL_OPTIONS="--no-print-usage" \
  -DCMAKE_CXX_FLAGS="$ALL" -DCMAKE_CXX_FLAGS_DEBUG="$DEBUG" \
  $PREFIX 
cmake --build build $PAR

if [ -n "$PREFIX" ]; then
  cmake --install build
fi
