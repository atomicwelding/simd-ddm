vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO jkriege2/TinyTIFF
    REF "${VERSION}"
    SHA512 ee672e80f369c4297abbf288f55f9655990e5efd71cfe409ef25d85b740304e8b2fde13ee8aed122f12067a46522b8b0153dc9d14bb1663f83b2ac7a8d355dc2
    HEAD_REF master
)


vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DTinyTIFF_BUILD_TESTS=OFF
        -DTinyTIFF_BUILD_STATIC_LIBS=ON
        -DTinyTIFF_BUILD_SHARED_LIBS=ON
)

vcpkg_cmake_install()

vcpkg_cmake_config_fixup(PACKAGE_NAME "TinyTIFFShared" CONFIG_PATH "lib/cmake/TinyTIFFShared")

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/share")

file(INSTALL "${SOURCE_PATH}/LICENSE" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}" RENAME copyright)
configure_file("${CMAKE_CURRENT_LIST_DIR}/usage" "${CURRENT_PACKAGES_DIR}/share/${PORT}/usage" COPYONLY)