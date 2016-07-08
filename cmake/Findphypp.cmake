if(NOT PHYPP_FOUND)
    # fine phy++ headers
    find_path(PHYPP_INCLUDE_DIR phypp.hpp
        HINTS ${PHYPP_ROOT_DIR} PATH_SUFFIXES include)

    file(GLOB PHYPP_HEADERS ${PHYPP_INCLUDE_DIR}/phypp/*.hpp)

    find_path(PHYPP_COMPILER_DIR cphy++
        HINTS ${PHYPP_ROOT_DIR} PATH_SUFFIXES bin)

    find_path(PHYPP_REFGEN_DIR phy++-refgen
        HINTS ${PHYPP_ROOT_DIR} PATH_SUFFIXES bin)

    set(PHYPP_COMPILER ${PHYPP_COMPILER_DIR}/cphy++)
    set(PHYPP_REFGEN ${PHYPP_REFGEN_DIR}/phy++-refgen)

    mark_as_advanced(PHYPP_INCLUDE_DIR PHYPP_HEADERS PHYPP_REFGEN PHYPP_COMPILER)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(PHYPP DEFAULT_MSG
        PHYPP_INCLUDE_DIR PHYPP_HEADERS PHYPP_REFGEN PHYPP_COMPILER)

    set(PHYPP_INCLUDE_DIRS ${PHYPP_INCLUDE_DIR})

    # configure compilers
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.3))
            message(STATUS "clang version >= 3.3 (${CMAKE_CXX_COMPILER_VERSION})")
        else()
            message(FATAL_ERROR "phy++ requires advanced features from the C++11 norm that are only available with clang 3.3 or higher (your version: ${CMAKE_CXX_COMPILER_VERSION}). Please upgrade your compiler.")
        endif()

        add_definitions(-Weverything)
        add_definitions(-Wno-c++98-compat-pedantic)
        add_definitions(-Wno-c++98-compat)
        add_definitions(-Wno-unused-parameter)
        add_definitions(-Wno-sign-conversion)
        add_definitions(-Wno-conversion)
        add_definitions(-Wno-missing-variable-declarations)
        add_definitions(-Wno-missing-prototypes)
        add_definitions(-Wno-padded)
        add_definitions(-Wno-float-equal)
        add_definitions(-Wno-unused-variable)
        add_definitions(-Wno-global-constructors)
        add_definitions(-Wno-exit-time-destructors)
        add_definitions(-Wno-weak-vtables)
        add_definitions(-Wno-covered-switch-default)
        add_definitions(-Wno-documentation-unknown-command)
        add_definitions(-Wno-unneeded-internal-declaration)
        add_definitions(-Wno-unused-function)
        add_definitions(-Wno-unused-macros)

        if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5))
            add_definitions(-Wno-old-style-cast)
        endif()

        add_definitions(-std=c++11)
        add_definitions(-ftemplate-backtrace-limit=0)
        add_definitions(-ferror-limit=5)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        if ("${CMAKE_CXX_COMPILER_VERSION}" STREQUAL "")
            message(WARNING "could not figure out the version of gcc, let's hope it is >= 4.7")
        else()
            if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7))
                message(STATUS "gcc version >= 4.7 (${CMAKE_CXX_COMPILER_VERSION})")
            else()
                message(FATAL_ERROR "phy++ requires advanced features from the C++11 norm that are only available with gcc 4.7 or higher (your version: ${CMAKE_CXX_COMPILER_VERSION}). Please upgrade your compiler.")
            endif()
        endif()

        add_definitions(-Wall)
        add_definitions(-std=c++11)
        add_definitions(-fmax-errors=5)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      message(ERROR "Intel C++ compiler is not supported")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
      message(ERROR "Microsoft Visual C++ compiler is not supported")
    endif()

    # find required libraries
    find_package(CFITSIO REQUIRED)
    find_package(Threads REQUIRED)

    set(PHYPP_INCLUDE_DIRS ${PHYPP_INCLUDE_DIRS} ${CFITSIO_INCLUDES})
    set(PHYPP_LIBRARIES ${PHYPP_LIBRARIES} ${CFITSIO_LIBRARIES})
    set(PHYPP_LIBRARIES ${PHYPP_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

    # find optional libraries
    if (NOT NO_LIBUNWIND)
        find_package(LibUnwind)
    endif()
    if (NOT NO_LIBDWARF)
        find_package(LibDwarf)
    endif()
    if (NOT NO_FFTW)
        find_package(FFTW 3)
    endif()
    if (NOT NO_PROFILER OR NOT NO_TCMALLOC)
        find_package(GooglePerfTools)
    endif()
    if (NOT NO_LAPACK)
        find_package(LAPACK)
    endif()
    if (NOT NO_GSL AND NOT NO_LAPACK)
        find_package(GSL)
    endif()
    if (NOT NO_WCSLIB)
        find_package(WCSLib)
    endif()

    # handle conditional reflection support
    # WIP: just disabled for now
    set(NO_REFLECTION 1)
    add_definitions(-DNO_REFLECTION)

    # handle conditional LAPACK support
    if (NOT LAPACK_FOUND OR NO_LAPACK)
        add_definitions(-DNO_LAPACK)
    else()
        set(PHYPP_LIBRARIES ${PHYPP_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()

    # handle conditional GSL support
    if (NOT GSL_FOUND OR NO_GSL OR NO_LAPACK OR NOT LAPACK_FOUND)
        add_definitions(-DNO_GSL)
    else()
        set(PHYPP_INCLUDE_DIRS ${PHYPP_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})
        set(PHYPP_LIBRARIES ${PHYPP_LIBRARIES} ${GSL_LIBRARIES})
    endif()

    # handle conditional WCSLib support
    if (NOT WCSLIB_FOUND OR NO_WCSLIB)
        add_definitions(-DNO_WCSLIB)
    else()
        set(PHYPP_INCLUDE_DIRS ${PHYPP_INCLUDE_DIRS} ${WCSLIB_INCLUDE_DIRS})
        set(PHYPP_LIBRARIES ${PHYPP_LIBRARIES} ${WCSLIB_LIBRARIES})
    endif()

    # handle conditional FFTW support
    if (NOT FFTW_FOUND OR NO_FFTW)
        add_definitions(-DNO_FFTW)
    else()
        set(PHYPP_INCLUDE_DIRS ${PHYPP_INCLUDE_DIRS} ${FFTW_INCLUDES})
        set(PHYPP_LIBRARIES ${PHYPP_LIBRARIES} ${FFTW_LIBRARIES})
    endif()

    # handle conditional LibUnwind support
    if (NOT LIBUNWIND_FOUND OR NO_LIBUNWIND)
        set(NO_UNWIND 1)
        add_definitions(-DNO_LIBUNWIND)
    else()
        set(NO_UNWIND 0)
        set(PHYPP_INCLUDE_DIRS ${PHYPP_INCLUDE_DIRS} ${LIBUNWIND_INCLUDE_DIR})
        set(PHYPP_LIBRARIES ${PHYPP_LIBRARIES} ${LIBUNWIND_LIBRARIES})
    endif()

    # handle conditional LibDwarf support
    if (NO_UNWIND OR NOT LIBDWARF_FOUND OR NO_LIBDWARF)
        add_definitions(-DNO_LIBDWARF)
    else()
        set(PHYPP_INCLUDE_DIRS ${PHYPP_INCLUDE_DIRS} ${LIBDWARF_INCLUDE_DIRS})
        set(PHYPP_LIBRARIES ${PHYPP_LIBRARIES} ${LIBDWARF_LIBRARIES})
    endif()

    # handle conditional Google perftools support
    if (TCMALLOC_LIBRARY)
        set(PHYPP_LIBRARIES ${PHYPP_LIBRARIES} ${TCMALLOC_LIBRARY})
    endif()
endif()

# Function to compile a phy++ program in CMake
# This feature does not support including/linking other external libraries, as well as
# "#define" commands, and is therefore not very powerful. But it is sufficient for basic
# needs. Supports reflection.
function(add_phypp_target CPP_FILE_NAME)
    # Generate binary name from c++ file
    get_filename_component(FILE_BASE ${CPP_FILE_NAME} NAME_WE)

    # Define the command to generate the binary file
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_custom_command(OUTPUT "${FILE_BASE}-make" VERBATIM COMMAND
            ${PHYPP_COMPILER} debug "${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME}" -o "${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make"
            DEPENDS ${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME} ${PHYPP_HEADERS})
    else()
        add_custom_command(OUTPUT "${FILE_BASE}-make" VERBATIM COMMAND
            ${PHYPP_COMPILER} optimize "${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME}" -o "${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make"
            DEPENDS ${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME} ${PHYPP_HEADERS})
    endif()

    # Create a target that will call this command
    add_custom_target(${FILE_BASE} ALL DEPENDS
        ${PROJECT_SOURCE_DIR}/${CPP_FILE_NAME}
        ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make
        ${PHYPP_HEADERS})

    # Specify installation of the binary
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}-make
        DESTINATION bin RENAME ${FILE_BASE})
endfunction()
