cmake_minimum_required(VERSION 2.6)
project(FAST++ C CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to RelWithDebInfo (-DCMAKE_BUILD_TYPE=...)")
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
find_package(phypp)

set(FASTPP_SHARE_DIR ${CMAKE_INSTALL_PREFIX}/share/FAST++)
install(FILES
    share/FILTER.RES.latest
    DESTINATION ${FASTPP_SHARE_DIR} COMPONENT data)

add_definitions(-DFASTPP_SHARE_DIR="${FASTPP_SHARE_DIR}")

include_directories(${PHYPP_INCLUDE_DIRS})

add_executable(fast++ src/fast++-read_input.cpp src/fast++-gridder.cpp src/fast++-fitter.cpp src/fast++.cpp)
target_link_libraries(fast++ ${PHYPP_LIBRARIES})
install(TARGETS fast++ DESTINATION bin)
