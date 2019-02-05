find_path(FMTLIB_INCLUDE_DIR fmt/format.h HINTS ${PROJECT_SOURCE_DIR}/External/include CMAKE_FIND_ROOT_PATH_BOTH)

find_package(PackageHandleStandardArgs)

find_package_handle_standard_args(fmtlib
    REQUIRED_VARS FMTLIB_INCLUDE_DIR
)

if(fmtlib_FOUND)
    set(FMTLIB_SOURCES
        ${FMTLIB_INCLUDE_DIR}/fmt/format.cc
        ${FMTLIB_INCLUDE_DIR}/fmt/posix.cc
    )
endif()
