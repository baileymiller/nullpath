# Module to find Field3D

# This module will first look into the directories defined by the variables:
#   - FIELD3D_HOME
#
# To use a custom OpenEXR
#   - Set the variable OPENEXR_CUSTOM to True
#   - Set the variable OPENEXR_CUSTOM_LIBRARY to the name of the library to
#     use, e.g. "SpiIlmImf"
#
# This module defines the following variables:
#
# FIELD3D_INCLUDE_DIR - where to find ImfRgbaFile.h, OpenEXRConfig, etc.
# FIELD3D_LIBRARY   - list of libraries to link against when using OpenEXR.
#                       This list does NOT include the IlmBase libraries.
#                       These are defined by the FindIlmBase module.
# FIELD3D_FOUND       - True if Field3D was found.

IF (WIN32)
	FIND_PATH( FIELD3D_INCLUDE_DIR Field3D/Field3DFile.h
		$ENV{FIELD3D_ROOT}/include/
		${FIELD3D_HOME}/include/
		DOC "The directory where Field3DFile.h resides")
	FIND_LIBRARY( FIELD3D_LIBRARY_RELEASE
		NAMES Field3D
		PATHS
		$ENV{FIELD3D_ROOT}/lib/
		${FIELD3D_HOME}/lib/    
		DOC "The OpenVDB library")
	FIND_LIBRARY( FIELD3D_LIBRARY_DEBUG
		NAMES Field3D_D
		PATHS
		$ENV{FIELD3D_ROOT}/lib/
		${FIELD3D_HOME}/lib    
		DOC "The OpenVDB library")
ELSE (WIN32)
        message("Searching for field3d include file in :" $ENV{FIELD3D_HOME})
        FIND_PATH(FIELD3D_INCLUDE_DIR Field3D/Field3DFile.h 
		PATHS $ENV{FIELD3D_HOME}/include/
                DOC "The directory where field3dfile.h resides"
                )
	FIND_LIBRARY(FIELD3D_LIBRARY_RELEASE NAMES libField3D.so Field3D libField3D 
		PATHS $ENV{FIELD3D_HOME}/lib/
                DOC "The field3d library")

	FIND_LIBRARY(FIELD3D_LIBRARY_DEBUG NAMES Field3D libField3D libField3D.so
		PATHS $ENV{FIELD3D_HOME}/lib/
                DOC "The field3d library")
ENDIF (WIN32)

UNSET ( FIELD3D_FOUND CACHE )

IF (NOT EXISTS ${FIELD3D_INCLUDE_DIR} )   
   SET ( FIELD3D_INCLUDE_DIR "" )
endif ()

IF (NOT EXISTS ${FIELD3D_LIBRARY_RELEASE} )   
   SET ( FIELD3D_LIBRARY_RELEASE "" )
endif ()

IF (NOT EXISTS ${FIELD3D_LIBRARY_DEBUG})
	SET (FIELD3D_LIBRARY_DEBUG "")
endif()

UNSET ( FIELD3D_FOUND CACHE)

# We just check for release libraries.
IF ( (NOT ${FIELD3D_INCLUDE_DIR} STREQUAL "") AND (NOT ${FIELD3D_LIBRARY_RELEASE} STREQUAL "") )
	SET( FIELD3D_FOUND TRUE CACHE BOOL "")
  message ( STATUS "Field3D Library: Found at ${FIELD3D_LIBRARY}" )
ELSE ()
  SET ( FIELD3D_FOUND FALSE CACHE BOOL "")
  message ( "Field3D Library: Not found. Confirm that FIELD3D_HOME is correct. Set FIELD3D_ROOT folder to correct directory.." )
ENDIF ()

MARK_AS_ADVANCED( FIELD3D_FOUND )
