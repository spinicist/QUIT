vcpkg_buildpath_length_warning(37)

block(SCOPE_FOR VARIABLES PROPAGATE SOURCE_PATH)
set(VCPKG_BUILD_TYPE release) # header-only

vcpkg_from_gitlab(
    GITLAB_URL https://gitlab.com
    OUT_SOURCE_PATH SOURCE_PATH
    REPO libeigen/eigen
    REF 14f845a1a87b77687cd328300d68733f158944ce
    SHA512 e1321109f5a0c250c47a320494e9835ba7448dcbdf3d8015e25cf99ec5f6449a0f2f3a7a1d2e51e0688ff57999d417a5a6f4a0721cdd24ddd54162a1aace7b84
    HEAD_REF master
    PATCHES
)

# vcpkg_from_git(
#     URL file:///Users/tobias/Code/eigen
#     OUT_SOURCE_PATH SOURCE_PATH
#     REF 78cb63a5d465d708be850c011d450d92d44d8434
# )

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DBUILD_TESTING=OFF
        -DEIGEN_BUILD_DOC=OFF
        -DEIGEN_BUILD_PKGCONFIG=ON
        -DEIGEN_BUILD_BLAS=OFF
        -DEIGEN_BUILD_LAPACK=OFF
        -DEIGEN_BUILD_DEMOS=OFF
        "-DCMAKEPACKAGE_INSTALL_DIR=${CURRENT_PACKAGES_DIR}/share/eigen3"
        "-DPKGCONFIG_INSTALL_DIR=${CURRENT_PACKAGES_DIR}/lib/pkgconfig"
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup()
vcpkg_replace_string("${CURRENT_PACKAGES_DIR}/share/eigen3/Eigen3Config.cmake" "if (NOT TARGET eigen)" "if (NOT TARGET Eigen3::Eigen)")
endblock()

if(NOT VCPKG_BUILD_TYPE)
    file(INSTALL "${CURRENT_PACKAGES_DIR}/lib/pkgconfig/eigen3.pc" DESTINATION "${CURRENT_PACKAGES_DIR}/debug/lib/pkgconfig")
endif()
vcpkg_fixup_pkgconfig()

file(GLOB INCLUDES "${CURRENT_PACKAGES_DIR}/include/eigen3/*")
# Copy the eigen header files to conventional location for user-wide MSBuild integration
file(COPY ${INCLUDES} DESTINATION "${CURRENT_PACKAGES_DIR}/include")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/COPYING.README")
