cmake_minimum_required(VERSION 2.6)
project(fast++ C CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release (-DCMAKE_BUILD_TYPE=...)")
    set(CMAKE_BUILD_TYPE "Release")
endif()

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
    fast++-read_input.cpp
    fast++-gridder.cpp
    fast++-gridder-ised.cpp
    fast++-gridder-custom.cpp
    fast++-fitter.cpp
    fast++-write_output.cpp
    fast++.cpp)
target_link_libraries(fast++ ${PHYPP_LIBRARIES})
target_link_libraries(fast++ ${TINYEXPR_LIBRARY})
install(TARGETS fast++ DESTINATION bin)

# Build FAST++ helper tools
add_executable(fast++-grid2fits fast++-grid2fits.cpp)
target_link_libraries(fast++-grid2fits ${PHYPP_LIBRARIES})
install(TARGETS fast++-grid2fits DESTINATION bin)
