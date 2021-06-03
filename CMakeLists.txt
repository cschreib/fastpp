cmake_minimum_required(VERSION 2.6)
project(FAST++-superbuild C CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected (-DCMAKE_BUILD_TYPE=...)")
    message(STATUS "Default to Release.")
    set(CMAKE_BUILD_TYPE "Release")
endif()

include(ExternalProject)

# Set default install dir to "here"
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(STATUS "No install directory selected (-DCMAKE_INSTALL_PREFIX=...)")
    message(STATUS "Default to ${PROJECT_SOURCE_DIR}/bin")
    set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" CACHE PATH "default install path" FORCE)
endif()

if (CMAKE_INSTALL_PREFIX STREQUAL ${PROJECT_SOURCE_DIR})
    set(FASTPP_SHARE_DIR ${PROJECT_SOURCE_DIR}/share)
    set(INSTALL_SHARED_DATA "0")
else()
    set(FASTPP_SHARE_DIR ${CMAKE_INSTALL_PREFIX}/share/FAST++)
    set(INSTALL_SHARED_DATA "1")
endif()

# Configure and build cfitsio
ExternalProject_Add(cfitsio
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/cfitsio
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external/cfitsio
    GIT_REPOSITORY https://github.com/cschreib/cfitsio.git
    GIT_TAG 9a9fbdba78f5583795f995c0f43007c191ff1b91
    CMAKE_ARGS
        -DBUILD_SHARED_LIBS=OFF
        -DUSE_PTHREADS=ON
        -DBUILD_TOOLS=OFF
        -DBUILD_TESTS=OFF
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_BUILD_TYPE=Release)

# Configure and build vif
ExternalProject_Add(vif
    DEPENDS cfitsio
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/vif
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external/vif
    GIT_REPOSITORY https://github.com/cschreib/vif.git
    GIT_TAG 0297ae9dde2dd13d983e1a86da189a95066fe7c7
    CMAKE_ARGS
        -DNO_REFLECTION=ON
        -DNO_FFTW=ON
        -DNO_LAPACK=ON
        -DNO_GSL=ON
        -DNO_WCSLIB=ON
        -DVIF_INPLACE_BUILD=ON
        -DCFITSIO_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_BUILD_TYPE=Release)

# Configure and build tinyexpr
ExternalProject_Add(tinyexpr
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/tinyexpr
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external/tinyexpr
    GIT_REPOSITORY https://github.com/cschreib/tinyexpr
    GIT_TAG d1ec9e328ba0be3743cce1856338f9a67fb608ef
    CMAKE_ARGS
        -DTE_NAT_LOG=ON
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_BUILD_TYPE=Release)

# Configure and build FAST++
ExternalProject_Add(fast++
    DEPENDS vif tinyexpr
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/src
    CMAKE_ARGS
        -DCFITSIO_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DVIF_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DTINYEXPR_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/external/include
        -DTINYEXPR_LIBRARY=${PROJECT_SOURCE_DIR}/external/lib/libtinyexpr.a
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}
        -DCMAKE_MODULE_PATH=${PROJECT_SOURCE_DIR}/external/vif/cmake
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DFASTPP_SHARE_DIR=${FASTPP_SHARE_DIR})

if (INSTALL_SHARED_DATA)
    install(FILES
        share/FILTER.RES.latest
        share/TEMPLATE_ERROR.fast.v0.2
        DESTINATION ${FASTPP_SHARE_DIR} COMPONENT data)
endif()

install(PROGRAMS
    ${CMAKE_BINARY_DIR}/bin/fast++
    ${CMAKE_BINARY_DIR}/bin/fast++-grid2fits
    ${CMAKE_BINARY_DIR}/bin/fast++-sfh2sed
    DESTINATION bin COMPONENT runtime)
