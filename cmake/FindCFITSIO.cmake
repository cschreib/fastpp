# +-----------------------------------------------------------------------------+
# | Copyright (C) 2011 |
# | Lars B"ahren (lbaehren@gmail.com) |
# | |
# | This program is free software; you can redistribute it and/or modify |
# | it under the terms of the GNU General Public License as published by |
# | the Free Software Foundation; either version 2 of the License, or |
# | (at your option) any later version. |
# | |
# | This program is distributed in the hope that it will be useful, |
# | but WITHOUT ANY WARRANTY; without even the implied warranty of |
# | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the |
# | GNU General Public License for more details. |
# | |
# | You should have received a copy of the GNU General Public License |
# | along with this program; if not, write to the |
# | Free Software Foundation, Inc., |
# | 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. |
# +-----------------------------------------------------------------------------+

# - Check for the presence of CFITSIO
#
# The following variables are set when CFITSIO is found:
# CFITSIO_FOUND = Set to true, if all components of CFITSIO have been found.
# CFITSIO_INCLUDES = Include path for the header files of CFITSIO
# CFITSIO_LIBRARIES = Link these to use CFITSIO
# CFITSIO_LFLAGS = Linker flags (optional)

if (NOT CFITSIO_FOUND)

  if (NOT CFITSIO_ROOT_DIR)
    set (CFITSIO_ROOT_DIR ${CMAKE_INSTALL_PREFIX})
  endif (NOT CFITSIO_ROOT_DIR)

  ##_____________________________________________________________________________
  ## Check for the header files

  find_path (CFITSIO_INCLUDES fitsio.h fitsio2.h
    HINTS ${CFITSIO_ROOT_DIR}
    PATHS /sw /usr /usr/local /opt/local
    PATH_SUFFIXES include include/fitsio include/cfitsio
    )

  ##_____________________________________________________________________________
  ## Check for the library

  find_library (CFITSIO_LIBRARIES cfitsio
    HINTS ${CFITSIO_ROOT_DIR}
    PATHS /sw /usr /usr/local /opt/local
    PATH_SUFFIXES lib
    )

  ##_____________________________________________________________________________
  ## Actions taken when all components have been found

  if (CFITSIO_INCLUDES AND CFITSIO_LIBRARIES)
    set (CFITSIO_FOUND TRUE)
  else (CFITSIO_INCLUDES AND CFITSIO_LIBRARIES)
    set (CFITSIO_FOUND FALSE)
    if (NOT CFITSIO_FIND_QUIETLY)
      if (NOT CFITSIO_INCLUDES)
        message (STATUS "Unable to find CFITSIO header files!")
      endif (NOT CFITSIO_INCLUDES)
      if (NOT CFITSIO_LIBRARIES)
        message (STATUS "Unable to find CFITSIO library files!")
      endif (NOT CFITSIO_LIBRARIES)
    endif (NOT CFITSIO_FIND_QUIETLY)
  endif (CFITSIO_INCLUDES AND CFITSIO_LIBRARIES)

  if (CFITSIO_FOUND)
    if (NOT CFITSIO_FIND_QUIETLY)
      message (STATUS "Found components for CFITSIO")
      message (STATUS "CFITSIO_INCLUDES = ${CFITSIO_INCLUDES}")
      message (STATUS "CFITSIO_LIBRARIES = ${CFITSIO_LIBRARIES}")
    endif (NOT CFITSIO_FIND_QUIETLY)
  else (CFITSIO_FOUND)
    if (CFITSIO_FIND_REQUIRED)
      message (FATAL_ERROR "Could not find CFITSIO!")
    endif (CFITSIO_FIND_REQUIRED)
  endif (CFITSIO_FOUND)

  ##_____________________________________________________________________________
  ## Mark advanced variables

  mark_as_advanced (
    CFITSIO_INCLUDES
    CFITSIO_LIBRARIES
    )

endif (NOT CFITSIO_FOUND)
