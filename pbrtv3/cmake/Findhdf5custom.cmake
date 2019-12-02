# Module to find HDF5 for field3d

# This module will first look into the directories defined by the variables:
#   - HDF5_ROOT
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
	FIND_PATH( HDF5_INCLUDE_DIR hdf5.h $ENV{HDF5_ROOT}/include/ ${HDF5_HOME}/include/ DOC "The directory where hdf5.h resides")
	FIND_LIBRARY( HDF5_LIBRARY_RELEASE NAMES libhdf5 PATHS $ENV{HDF5_ROOT}/lib/ ${HDF5_HOME}/lib DOC "The hdf5 library")
	FIND_LIBRARY( HDF5_CPP_LIBRARY_RELEASE NAMES libhdf5_cpp PATHS $ENV{HDF5_ROOT}/lib/ ${HDF5_HOME}/lib DOC "The hdf5_cpp library")
	FIND_LIBRARY( SZIP_LIBRARY_RELEASE NAMES libszip PATHS $ENV{HDF5_ROOT}/lib/ ${HDF5_HOME}/lib DOC "The hdf5 szip library")
	FIND_LIBRARY( ZLIB_LIBRARY_RELEASE NAMES libzlib PATHS $ENV{HDF5_ROOT}/lib/ ${HDF5_HOME}/lib DOC "The hdf5 zlib library")
	
	FIND_LIBRARY( HDF5_LIBRARY_DEBUG NAMES libhdf5_d PATHS $ENV{HDF5_ROOT}/lib/ ${HDF5_HOME}/lib DOC "The hdf5 debug library")
	FIND_LIBRARY( HDF5_CPP_LIBRARY_DEBUG NAMES libhdf5_cpp_d PATHS $ENV{HDF5_ROOT}/lib/ ${HDF5_HOME}/lib DOC "The hdf5_cpp debug library")
	FIND_LIBRARY( SZIP_LIBRARY_DEBUG NAMES libszip_d PATHS $ENV{HDF5_ROOT}/lib/ ${HDF5_HOME}/lib DOC "The hdf5 szip debug library")
	FIND_LIBRARY( ZLIB_LIBRARY_DEBUG NAMES libzlib_d PATHS $ENV{HDF5_ROOT}/lib/ ${HDF5_HOME}/lib DOC "The hdf5 zlib debug library")
ELSEIF (APPLE)
	message("findhdf5 custom searching within : " $ENV{HDF5_HOME})
	FIND_PATH( HDF5_INCLUDE_DIR hdf5.h {HDF5_ROOT}/include/ $ENV{HDF5_HOME}/include/ DOC "The directory where hdf5.h resides")
	FIND_LIBRARY( HDF5_LIBRARY_RELEASE NAMES libhdf5.dylib hdf5 PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5 library")
	FIND_LIBRARY( HDF5_CPP_LIBRARY_RELEASE NAMES libhdf5_cpp.dylib hdf5_cpp PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5_cpp library")
	FIND_LIBRARY( SZIP_LIBRARY_RELEASE NAMES libsz.dylib sz PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5 szip library")
	FIND_LIBRARY( ZLIB_LIBRARY_RELEASE NAMES libz.dylib z PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5 zlib library")
	
	FIND_LIBRARY( HDF5_LIBRARY_DEBUG NAMES libhdf5.dylib hdf5 PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5 debug library")
	FIND_LIBRARY( HDF5_CPP_LIBRARY_DEBUG NAMES libhdf5_cpp.dylib hdf5_cpp PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5_cpp debug library")
	FIND_LIBRARY( SZIP_LIBRARY_DEBUG NAMES libsz.dylib sz PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5 szip debug library")
	FIND_LIBRARY( ZLIB_LIBRARY_DEBUG NAMES libz.dylib z PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5 zlib debug library")
ELSE ()
	message("findhdf5 custom searching within : " $ENV{HDF5_HOME})
	FIND_PATH( HDF5_INCLUDE_DIR hdf5.h 
		{HDF5_ROOT}/include/hdf5/serial 
		$ENV{HDF5_HOME}/include/hdf5/serial 
		$ENV{HDF5_HOME}/include/
		DOC "The directory where hdf5.h resides")
	message("HDF5 include dir : " ${HDF5_INCLUDE_DIR})
	FIND_LIBRARY( HDF5_LIBRARY_RELEASE NAMES libhdf5_serial.dylib hdf5_serial hdf5 PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5 library")
	FIND_LIBRARY( HDF5_CPP_LIBRARY_RELEASE NAMES libhdf5_cpp.dylib hdf5_cpp hdf5_cpp PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5_cpp library")
	FIND_LIBRARY( SZIP_LIBRARY_RELEASE NAMES libsz.dylib sz PATHS $ENV{HD5_ROOT}/lib/ $ENV{HDF5_HOME}/lib $ENV{SZIP_HOME}/lib/ DOC "The hdf5 szip library")
	FIND_LIBRARY( ZLIB_LIBRARY_RELEASE NAMES libz.dylib z PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5 zlib library")
	
	FIND_LIBRARY( HDF5_LIBRARY_DEBUG NAMES libhdf5_serial.dylib hdf5_serial hdf5 PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5 debug library")
	FIND_LIBRARY( HDF5_CPP_LIBRARY_DEBUG NAMES libhdf5_cpp.dylib hdf5_cpp hdf5_cpp PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5_cpp debug library")
	FIND_LIBRARY( SZIP_LIBRARY_DEBUG NAMES libsz.dylib sz PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib $ENV{SZIP_HOME}/lib/ DOC "The hdf5 szip debug library")
	FIND_LIBRARY( ZLIB_LIBRARY_DEBUG NAMES libz.dylib z PATHS $ENV{HDF5_ROOT}/lib/ $ENV{HDF5_HOME}/lib DOC "The hdf5 zlib debug library")
ENDIF ()

UNSET ( HDF5_FOUND CACHE )

IF (NOT EXISTS ${HDF5_INCLUDE_DIR} )   
   SET ( HDF5_INCLUDE_DIR "" )
endif ()

IF (NOT EXISTS ${HDF5_LIBRARY_RELEASE} )   
   SET ( HDF5_LIBRARY_RELEASE "" )
endif ()

IF (NOT EXISTS ${HDF5_LIBRARY_DEBUG})
	SET (HDF5_LIBRARIES_DEBUG "")
endif()

UNSET ( HDF5_FOUND CACHE)

# We just check for release libraries.
IF ( (NOT ${HDF5_INCLUDE_DIR} STREQUAL "") AND (NOT ${HDF5_LIBRARY_RELEASE} STREQUAL "") )
	SET( HDF5_FOUND TRUE CACHE BOOL "")
  message ( STATUS "hdf5 Library: Found at ${HDF5_LIBRARY_RELEASE}" )
ELSE ()
  SET ( HDF5_FOUND FALSE CACHE BOOL "")
  message ( "hdf5 Library: Not found. Confirm that hdf5_ROOT is correct. Set HDF5_ROOT folder to correct directory.." )
ENDIF ()

MARK_AS_ADVANCED( HDF5_FOUND )
