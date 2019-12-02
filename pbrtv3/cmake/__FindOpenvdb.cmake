# - Find OPENVDB library
# Find the native OPENVDB includes and library
# This module defines
#  OPENVDB_INCLUDE_DIRS, where to find openvdb.h, Set when
#                            OPENVDB_INCLUDE_DIR is found.
#  OPENVDB_LIBRARIES, libraries to link against to use OPENVDB.
#  OPENVDB_ROOT_DIR, The base directory to search for OPENVDB.
#                        This can also be an environment variable.
#  OPENVDB_FOUND, If false, do not try to use OPENVDB.
#
# also defined, but not for general use are
#  OPENVDB_LIBRARY, where to find the OPENVDB library.

#=============================================================================
# Copyright 2015 Blender Foundation.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

# If OPENVDB_ROOT_DIR was defined in the environment, use it.
FIND_PATH(OPENVDB_INCLUDE_DIR openvdb/openvdb.h $ENV{OPENVDB_ROOT}/include/ ${OPENVDB_HOME}/include/ DOC "The directory where openvdb.h resides")
FIND_LIBRARY(OPENVDB_LIBRARY_DEBUG NAMES openvdb PATHS $ENV{FIELD3D_ROOT}/lib/ ${FIELD3D_HOME}/lib DOC "The OpenVDB debug static library ")
FIND_LIBRARY(OPENVDB_LIBRARY_RELEASE NAMES openvdb PATHS $ENV{FIELD3D_ROOT}/lib/ ${FIELD3D_HOME}/lib DOC "The OpenVDB release static library ")


# handle the QUIETLY and REQUIRED arguments and set OPENVDB_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OPENVDB DEFAULT_MSG
    OPENVDB_LIBRARY_DEBUG OPENVDB_LIBRARY_RELEASE OPENVDB_INCLUDE_DIR)

IF(OPENVDB_FOUND)
  SET(OPENVDB_LIBRARIES_DEBUG ${OPENVDB_LIBRARY_DEBUG})
  SET(OPENVDB_LIBRARIES_RELEASE ${OPENVDB_LIBRARY_RELEASE})
  SET(OPENVDB_INCLUDE_DIRS ${OPENVDB_INCLUDE_DIR})
ENDIF(OPENVDB_FOUND)

MARK_AS_ADVANCED(
  OPENVDB_INCLUDE_DIR
  OPENVDB_LIBRARY
)

UNSET(_openvdb_SEARCH_DIRS)
