# - Try to find LAPACK
# Variables used by this module:
#  LAPACK_ROOT_DIR     - LAPACK root directory
# Variables defined by this module:
#  LAPACK_FOUND        - system has LAPACK
#  LAPACK_LIBRARY      - the LAPACK library (cached)
#  LAPACK_LIBRARIES    - the LAPACK libraries

if(NOT LAPACK_FOUND)

  find_library(LAPACK_LIBRARY lapack
    HINTS ${LAPACK_ROOT_DIR} PATH_SUFFIXES lib)
  mark_as_advanced(LAPACK_INCLUDE_DIR LAPACK_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(LAPACK DEFAULT_MSG
    LAPACK_LIBRARY)

  set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})

endif(NOT LAPACK_FOUND)
