# Module to find OpenVDB

# This module will first look into the directories defined by the variables:
#   - OPENVDB_HOME
#
# To use a custom OpenEXR
#   - Set the variable OPENEXR_CUSTOM to True
#   - Set the variable OPENEXR_CUSTOM_LIBRARY to the name of the library to
#     use, e.g. "SpiIlmImf"
#
# This module defines the following variables:
#
# OPENVDB_INCLUDE_DIR - where to find ImfRgbaFile.h, OpenEXRConfig, etc.
# OPENVDB_LIBRARIES   - list of libraries to link against when using OpenEXR.
#                       This list does NOT include the IlmBase libraries.
#                       These are defined by the FindIlmBase module.
# OPENVDB_FOUND       - True if OpenEXR was found.
IF (WIN32)
	FIND_PATH( OPENVDB_INCLUDE_DIR openvdb.h
		$ENV{OPENVDB_ROOT}/include/openvdb/
		${OPENVDB_HOME}/include/
		DOC "The directory where openvdb.h resides")
	FIND_LIBRARY( OPENVDB_LIBRARY
		NAMES openvdb
		PATHS
		$ENV{OPENVDB_ROOT}/lib/
		${OPENVDB_HOME}/lib    
		DOC "The OpenVDB library")
ELSE ()
	FIND_PATH(OPENVDB_INCLUDE_DIR openvdb.h
                $ENV{OPENVDB_ROOT}/include/openvdb/
		${OPENVDB_HOME}/include/
		/usr/include
		/usr/local/include
		/sw/include
		/opt/local/include
		DOC "The directory where openvdb.h resides")
	FIND_LIBRARY( OPENVDB_LIBRARY
		NAMES openvdb
		PATHS
		/usr/lib64
		/usr/lib
		/usr/local/lib64
		/usr/local/lib
		/sw/lib
		/opt/local/lib
                $ENV{OPENVDB_ROOT}/lib/
		${OPENVDB_HOME}/lib    
		DOC "The OpenVDB library")
ENDIF ()

UNSET ( OPENVDB_FOUND CACHE )

IF (NOT EXISTS ${OPENVDB_INCLUDE_DIR} )   
   SET ( OPENVDB_INCLUDE_DIR "" )
endif ()

IF (NOT EXISTS ${OPENVDB_LIBRARY} )   
   SET ( OPENVDB_LIBRARY "" )
endif ()

UNSET ( OPENVDB_FOUND CACHE)

IF ( (NOT ${OPENVDB_INCLUDE_DIR} STREQUAL "") AND (NOT ${OPENVDB_LIBRARY} STREQUAL "") )
	SET( OPENVDB_FOUND TRUE CACHE BOOL "")
  message ( STATUS "OpenVDB Library: Found at ${OPENVDB_LIBRARY}" )
ELSE ()
  SET ( OPENVDB_FOUND FALSE CACHE BOOL "")
  message ( "OpenVDB Library: Not found. Confirm that OPENVDB_HOME is correct." )
ENDIF ()

MARK_AS_ADVANCED( OPENVDB_FOUND )
