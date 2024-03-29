find_package(Git REQUIRED)

if(WIN32)
    include(InstallRequiredSystemLibraries)
    set(CPACK_GENERATOR NSIS)
else()
    set(CPACK_GENERATOR TGZ)
endif()

set(CPACK_PACKAGE_NAME cryptoutils)
set(CPACK_PACKAGE_VENDOR arxflay)
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/LICENSE)
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CPACK_PACKAGE_VERSION_MAJOR ${VER_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VER_MINOR})
set(CPACK_PACKAGE_VERSION ${VER_MAJOR}.${VER_MINOR})
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${PKG_VER_NAME}")
set(CPACK_PACKAGE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/installers/${CPACK_GENERATOR}")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Alexej Fedorenko(arxflay)")
#this must be a last entry in this file
include(CPack)
