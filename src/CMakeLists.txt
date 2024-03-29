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

find_package(vif REQUIRED)

# Figure out git hash, if any
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  RESULT_VARIABLE GIT_COMMAND_SUCCESS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if (GIT_COMMAND_SUCCESS EQUAL 0)
    add_definitions(-DFASTPP_GIT_HASH="${GIT_COMMIT_HASH}")
else()
    add_definitions(-DFASTPP_GIT_HASH="")
endif()

# Setup share directory
add_definitions(-DFASTPP_SHARE_DIR="${FASTPP_SHARE_DIR}")

# Include dependencies
message(STATUS ${VIF_INCLUDE_DIRS})
include_directories(${VIF_INCLUDE_DIRS})
include_directories(${TINYEXPR_INCLUDE_DIR})

# Build FAST++
add_executable(fast++
    fast++-read_input.cpp
    fast++-ssp.cpp
    fast++-gridder.cpp
    fast++-gridder-ised.cpp
    fast++-gridder-custom.cpp
    fast++-fitter.cpp
    fast++-write_output.cpp
    fast++.cpp)
target_link_libraries(fast++ ${VIF_LIBRARIES})
target_link_libraries(fast++ ${TINYEXPR_LIBRARY})
install(TARGETS fast++ DESTINATION bin)

# Build FAST++ helper tools
add_executable(fast++-grid2fits fast++-grid2fits.cpp)
target_link_libraries(fast++-grid2fits ${VIF_LIBRARIES})
install(TARGETS fast++-grid2fits DESTINATION bin)

add_executable(fast++-sfh2sed fast++-sfh2sed.cpp fast++-ssp.cpp)
target_link_libraries(fast++-sfh2sed ${VIF_LIBRARIES})
install(TARGETS fast++-sfh2sed DESTINATION bin)
