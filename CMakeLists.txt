cmake_minimum_required(VERSION 2.6)
project(FAST++ C CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release (-DCMAKE_BUILD_TYPE=...)")
    set(CMAKE_BUILD_TYPE "Release")
endif()

include(ExternalProject)

# Configure and build cfitsio
ExternalProject_Add(cfitsio
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/cfitsio
    DOWNLOAD_COMMAND "" UPDATE_COMMAND ""
    CMAKE_ARGS
        -DBUILD_SHARED_LIBS=OFF
        -DUSE_PTHREADS=ON
        -DBUILD_TOOLS=OFF
        -DBUILD_TESTS=OFF
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_BUILD_TYPE=Release)

# Configure and build phy++
ExternalProject_Add(phy++
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/phypp
    DOWNLOAD_COMMAND "" UPDATE_COMMAND ""
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

set(CFITSIO_ROOT_DIR ${PROJECT_SOURCE_DIR}/external)
set(PHYPP_ROOT_DIR ${PROJECT_SOURCE_DIR}/external)
set(NO_REFLECTION ON)
set(NO_FFTW ON)
set(NO_LAPACK ON)
set(NO_GSL ON)
set(NO_WCSLIB ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/external/phypp/cmake")
find_package(phypp REQUIRED)

# Setup share directory
set(FASTPP_SHARE_DIR ${CMAKE_INSTALL_PREFIX}/share/FAST++)
install(FILES
    share/FILTER.RES.latest
    share/TEMPLATE_ERROR.fast.v0.2
    DESTINATION ${FASTPP_SHARE_DIR} COMPONENT data)

add_definitions(-DFASTPP_SHARE_DIR="${FASTPP_SHARE_DIR}")

# Include dependencies
include_directories(${PHYPP_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/external/include)

# Build FAST++
add_executable(fast++
    src/fast++-read_input.cpp
    src/fast++-gridder.cpp
    src/fast++-fitter.cpp
    src/fast++-write_output.cpp
    src/fast++.cpp)
target_link_libraries(fast++ ${PHYPP_LIBRARIES})
install(TARGETS fast++ DESTINATION bin)

add_dependencies(fast++ phy++)
add_dependencies(phy++ cfitsio)
