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
    GIT_REPOSITORY https://github.com/HEASARC/cfitsio.git
    UPDATE_COMMAND "" # No need to pull changes all the time, this is stable
    CMAKE_ARGS
        -DBUILD_SHARED_LIBS=OFF
        -DUSE_PTHREADS=ON
        -DBUILD_TOOLS=OFF
        -DBUILD_TESTS=OFF
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_BUILD_TYPE=Release)

# Configure and build phy++
ExternalProject_Add(phy++
    DEPENDS cfitsio
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/phypp
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external/phypp
    GIT_REPOSITORY https://github.com/cschreib/phypp.git
    UPDATE_COMMAND "" # No need to pull changes all the time, this is stable
    CMAKE_ARGS
        -DNO_REFLECTION=ON
        -DNO_FFTW=ON
        -DNO_LAPACK=ON
        -DNO_GSL=ON
        -DNO_WCSLIB=ON
        -DPHYPP_INPLACE_BUILD=ON
        -DCFITSIO_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_BUILD_TYPE=Release)

# Configure and build tinyexpr
ExternalProject_Add(tinyexpr
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/tinyexpr
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external/tinyexpr
    GIT_REPOSITORY https://github.com/cschreib/tinyexpr
    UPDATE_COMMAND "" # No need to pull changes all the time, this is stable
    CMAKE_ARGS
        -DTE_NAT_LOG=ON
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_BUILD_TYPE=Release)

# Point to external dependencies
set(CFITSIO_ROOT_DIR ${PROJECT_SOURCE_DIR}/external)
set(PHYPP_ROOT_DIR ${PROJECT_SOURCE_DIR}/external)
set(TINYEXPR_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external/include)
set(TINYEXPR_LIBRARY ${PROJECT_SOURCE_DIR}/external/lib/libtinyexpr.a)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/external/phypp/cmake)

# Setup phy++ options
set(NO_REFLECTION ON)
set(NO_FFTW ON)
set(NO_LAPACK ON)
set(NO_GSL ON)
set(NO_WCSLIB ON)

find_package(phypp REQUIRED)

# Setup share directory
add_definitions(-DFASTPP_SHARE_DIR="${FASTPP_SHARE_DIR}")

# Include dependencies
message(STATUS ${PHYPP_INCLUDE_DIRS})
include_directories(${PHYPP_INCLUDE_DIRS})
include_directories(${TINYEXPR_INCLUDE_DIR})

# Build FAST++
add_executable(fast++
    src/fast++-read_input.cpp
    src/fast++-ssp.cpp
    src/fast++-gridder.cpp
    src/fast++-gridder-ised.cpp
    src/fast++-gridder-custom.cpp
    src/fast++-fitter.cpp
    src/fast++-write_output.cpp
    src/fast++.cpp)
target_link_libraries(fast++ ${PHYPP_LIBRARIES})
target_link_libraries(fast++ ${TINYEXPR_LIBRARY})
install(TARGETS fast++ DESTINATION bin)

# Build FAST++ helper tools
add_executable(fast++-grid2fits src/fast++-grid2fits.cpp)
target_link_libraries(fast++-grid2fits ${PHYPP_LIBRARIES})
install(TARGETS fast++-grid2fits DESTINATION bin)

add_executable(fast++-sfh2sed src/fast++-sfh2sed.cpp src/fast++-ssp.cpp)
target_link_libraries(fast++-sfh2sed ${PHYPP_LIBRARIES})
install(TARGETS fast++-sfh2sed DESTINATION bin)

if (INSTALL_SHARED_DATA)
    install(FILES
        share/FILTER.RES.latest
        share/TEMPLATE_ERROR.fast.v0.2
        DESTINATION ${FASTPP_SHARE_DIR} COMPONENT data)
endif()
