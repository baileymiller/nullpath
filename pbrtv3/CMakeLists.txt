CMAKE_MINIMUM_REQUIRED ( VERSION 2.8 )
set(CMAKE_CXX_STANDARD 14)
# For sanitizers
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

PROJECT(PBRT-V3)

ENABLE_TESTING()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/ext/openexr/OpenEXR")
  message(FATAL_ERROR "The OpenEXR submodule directory is missing! "
    "You probably did not clone the project with --recursive. It is possible to recover "
    "by running \"git submodule update --init --recursive\"")
endif()

if(NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/ext/glog")
  message(FATAL_ERROR "The glog submodule directory is missing! "
    "You probably did not clone the project with --recursive, or first checked out "
    "pbrt before it was added. It is possible to recover by running "
    "\"git submodule update --init --recursive\"")
endif()

IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-conversion-null")
ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-register")
ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

  FIND_PROGRAM(XIAR xiar)
  IF(XIAR)
    SET(CMAKE_AR "${XIAR}")
  ENDIF(XIAR)
  MARK_AS_ADVANCED(XIAR)

  FIND_PROGRAM(XILD xild)
  IF(XILD)
    SET(CMAKE_LINKER "${XILD}")
  ENDIF(XILD)
  MARK_AS_ADVANCED(XILD)

  # ICC will default to -fp-model fast=1, which performs value-unsafe optimizations which will
  # cause pbrt_test to fail. For safety, -fp-model precise is explicitly set here by default.
  set(FP_MODEL "precise" CACHE STRING "The floating point model to compile with.")
  set_property(CACHE FP_MODEL PROPERTY STRINGS "precise" "fast=1" "fast=2")

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model ${FP_MODEL}")
ENDIF()

IF(MSVC)
  ADD_DEFINITIONS (/D _CRT_SECURE_NO_WARNINGS)
  ADD_DEFINITIONS (/D YY_NO_UNISTD_H)
ENDIF()

FIND_PACKAGE ( Sanitizers )
FIND_PACKAGE ( Threads )

# Add support for building with openmp for parallelization in some parts of the code
#FIND_PACKAGE(OpenMP REQUIRED)
#IF(OPENMP_FOUND)
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#ENDIF()


IF(CMAKE_BUILD_TYPE MATCHES RELEASE)
  ADD_DEFINITIONS (-DNDEBUG)
ENDIF()

# Optionally use Bison and Flex to regenerate parser files
# Use pregenerated files otherwise (may be outdated)
FIND_PACKAGE ( BISON )
FIND_PACKAGE ( FLEX )
IF(BISON_FOUND AND FLEX_FOUND)
  SET(BisonOutput ${CMAKE_BINARY_DIR}/pbrtparse.cpp)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${BisonOutput} ${CMAKE_BINARY_DIR}/pbrtparse.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/core/pbrtparse.y
    COMMAND ${BISON_EXECUTABLE} -d -v -t
            --output=${CMAKE_BINARY_DIR}/pbrtparse.c
            ${CMAKE_SOURCE_DIR}/src/core/pbrtparse.y
    COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_BINARY_DIR}/pbrtparse.c ${CMAKE_BINARY_DIR}/pbrtparse.cpp
    COMMENT "Generating pbrtparse.cpp and pbrtparse.h"
  )

  SET(FlexOutput ${CMAKE_BINARY_DIR}/pbrtlex.cpp)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${FlexOutput}
    DEPENDS ${CMAKE_BINARY_DIR}/pbrtparse.h ${CMAKE_SOURCE_DIR}/src/core/pbrtlex.ll
    COMMAND ${FLEX_EXECUTABLE}
            -o${FlexOutput}
            ${CMAKE_SOURCE_DIR}/src/core/pbrtlex.ll
    COMMENT "Generating pbrtlex.cpp"
  )

  SET ( PBRT_YACC_LEX_SOURCE
    ${BisonOutput}
    ${FlexOutput}
    )
ELSE()
  SET ( PBRT_YACC_LEX_SOURCE
    src/core/pbrtlex.cpp
    src/core/pbrtparse.cpp
    )
ENDIF()

IF(WIN32)
  # Build zlib (only on Windows)
  SET(ZLIB_BUILD_STATIC_LIBS ON CACHE BOOL " " FORCE)
  SET(ZLIB_BUILD_SHARED_LIBS OFF CACHE BOOL " " FORCE)
  ADD_SUBDIRECTORY(src/ext/zlib)

  SET(ZLIB_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/ext/zlib")
  SET(ZLIB_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/src/ext/zlib/$<CONFIGURATION>/zlibstatic.lib")
  SET_PROPERTY(TARGET zlibstatic PROPERTY FOLDER "dependencies")
  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR} "${CMAKE_CURRENT_BINARY_DIR}/src/ext/zlib")
ENDIF()

#SET(ILMBASE_NAMESPACE_VERSIONING OFF CACHE BOOL " " FORCE)
#SET(OPENEXR_NAMESPACE_VERSIONING OFF CACHE BOOL " " FORCE)
#SET(OPENEXR_BUILD_SHARED_LIBS    OFF CACHE BOOL " " FORCE)
#SET(ILMBASE_BUILD_SHARED_LIBS    OFF CACHE BOOL " " FORCE)
#ADD_SUBDIRECTORY(src/ext/openexr)
#SET_PROPERTY(TARGET IexMath eLut toFloat b44ExpLogTable dwaLookups CopyIlmBaseLibs IlmThread Half Iex Imath IlmImf PROPERTY FOLDER "dependencies")
#INCLUDE_DIRECTORIES (
#  src/ext/openexr/IlmBase/Imath
#  src/ext/openexr/IlmBase/Half
#  src/ext/openexr/IlmBase/Iex
#  src/ext/openexr/OpenEXR/IlmImf
#  ${CMAKE_BINARY_DIR}/src/ext/openexr/IlmBase/config
#  ${CMAKE_BINARY_DIR}/src/ext/openexr/OpenEXR/config
#)

IF(WIN32)
  SET(OPENEXR_LIBS IlmImf Imath Half zlibstatic)
ENDIF()

SET(WITH_GFLAGS OFF CACHE BOOL "Use gflags")
SET(BUILD_SHARED_LIBS OFF CACHE BOOL " " FORCE)
ADD_SUBDIRECTORY(src/ext/glog)
SET_PROPERTY(TARGET glog PROPERTY FOLDER "dependencies")
INCLUDE_DIRECTORIES (
  src/ext/glog/src
  ${CMAKE_BINARY_DIR}/src/ext/glog
)


# set the boost library dir variable if win32
#set(BOOST_ROOT /usr/local/Cellar/boost/HEAD-8111738_1/)
#set(Boost_INCLUDE_DIR /usr/local/Cellar/boost/HEAD-8111738_1/include/)
#set(BOOST_LIBRARYDIR /usr/local/Cellar/boost/HEAD-8111738_1/lib/)
#set(Boost_DIR /usr/local/Cellar/boost/HEAD-8111738_1/)
#set(Boost_USE_STATIC_LIBS   ON)

# Although pbrt comes with a openexr version, the paths are messed update
# hence we might needdylib)
FIND_PACKAGE(Field3D REQUIRED)
FIND_PACKAGE(OpenVDB REQUIRED)
FIND_PACKAGE(Boost REQUIRED COMPONENTS regex thread system iostreams)
FIND_PACKAGE(hdf5custom REQUIRED)
FIND_PACKAGE(TBB REQUIRED)
if (APPLE)
set(OPENEXR_LIBS 
      /usr/local/Cellar/openexr/2.2.0/lib/libIlmImf.dylib
      /usr/local/Cellar/ilmbase/2.2.0/lib/libHalf.dylib 
      /usr/local/Cellar/ilmbase/2.2.0/lib/libImath.dylib 
      /usr/local/Cellar/ilmbase/2.2.0/lib/libIex.dylib 
)
set(OPENEXR_INCLUDE_DIRS 
      /usr/local/Cellar/ilmbase/2.2.0/include/OpenEXR
      /usr/local/Cellar/openexr/2.2.0/include/OpenEXR
      /usr/local/Cellar/openexr/2.2.0/include/
)
elseif(UNIX)
FIND_PACKAGE(OpenEXR REQUIRED)
FIND_LIBRARY(ILMIMF_LIB 
		NAMES 	
			libIlmImf 
			IlmImf 
			libIlmImf.so
		PATHS 
			$ENV{OPENEXR_ROOT} 
			$ENV{ILMBASE_LIB_ROOT}/lib/$ENV{X86_LIB_EXT})
set(OPENEXR_LIBS 
	$ENV{ILMBASE_LIB_ROOT}/lib/$ENV{X86_LIB_EXT}/libHalf.so
	${ILMIMF_LIB}
	$ENV{ILMBASE_LIB_ROOT}/lib/$ENV{X86_LIB_EXT}/libImath.so
	$ENV{ILMBASE_LIB_ROOT}/lib/$ENV{X86_LIB_EXT}/libIex.so
	${OPENEXR_LIBS}
)

set(OPENEXR_INCLUDE_DIRS 
	$ENV{ILMBASE_INCLUDE_ROOT}/include/OpenEXR
	$ENV{OPENEXR_ROOT}/include/OpenEXR
	${OPENEXR_INCLUDE_DIRS}
)
endif()

# add include directories
include_directories(
	${Boost_INCLUDE_DIR}
  ${FIELD3D_INCLUDE_DIR}
	${HDF5_INCLUDE_DIR}
	${OPENEXR_INCLUDE_DIRS}
	${OPENVDB_INCLUDE_DIR}
	${TBB_INCLUDE_DIR}
  ${ILMBASE_INCLUDE_DIRS}
)



SET ( PBRT_CORE_SOURCE
  src/core/api.cpp
  src/core/bssrdf.cpp
  src/core/camera.cpp
  src/core/efloat.cpp
  src/core/error.cpp
  src/core/fileutil.cpp
  src/core/film.cpp
  src/core/filter.cpp
  src/core/floatfile.cpp
  src/core/geometry.cpp
  src/core/imageio.cpp
  src/core/integrator.cpp
  src/core/interaction.cpp
  src/core/interpolation.cpp
  src/core/light.cpp
  src/core/lightdistrib.cpp
  src/core/lowdiscrepancy.cpp
  src/core/material.cpp
  src/core/medium.cpp
  src/core/memory.cpp
  src/core/microfacet.cpp
  src/core/parallel.cpp
  src/core/paramset.cpp
  src/core/parser.cpp
  src/core/primitive.cpp
  src/core/progressreporter.cpp
  src/core/quaternion.cpp
  src/core/reflection.cpp
  src/core/sampler.cpp
  src/core/sampling.cpp
  src/core/scene.cpp
  src/core/shape.cpp
  src/core/sobolmatrices.cpp
  src/core/spectrum.cpp
  src/core/stats.cpp
  src/core/texture.cpp
  src/core/transform.cpp
  )

SET ( PBRT_CORE_HEADERS
  src/core/api.h
  src/core/bssrdf.h
  src/core/camera.h
  src/core/efloat.h
  src/core/error.h
  src/core/fileutil.h
  src/core/film.h
  src/core/filter.h
  src/core/floatfile.h
  src/core/geometry.h
  src/core/imageio.h
  src/core/integrator.h
  src/core/interaction.h
  src/core/interpolation.h
  src/core/light.h
  src/core/lowdiscrepancy.h
  src/core/material.h
  src/core/medium.h
  src/core/memory.h
  src/core/microfacet.h
  src/core/mipmap.h
  src/core/parallel.h
  src/core/paramset.h
  src/core/parser.h
  src/core/pbrt.h
  src/core/port.h
  src/core/primitive.h
  src/core/progressreporter.h
  src/core/quaternion.h
  src/core/reflection.h
  src/core/rng.h
  src/core/sampler.h
  src/core/sampling.h
  src/core/scene.h
  src/core/shape.h
  src/core/sobolmatrices.h
  src/core/spectrum.h
  src/core/stats.h
  src/core/stringprint.h
  src/core/texture.h
  src/core/transform.h
  )

FILE ( GLOB PBRT_SOURCE
  src/ext/*.c
  src/ext/*.cpp
  src/accelerators/*.cpp
  src/cameras/*.cpp
  src/filters/*.cpp
  src/integrators/*.cpp
  src/lights/*.cpp
  src/lights/sunsky/*.cpp
  src/materials/*.cpp
  src/samplers/*.cpp
  src/shapes/*.cpp
  src/textures/*.cpp
  src/media/*.cpp
  )

INCLUDE_DIRECTORIES ( src )
INCLUDE_DIRECTORIES ( src/core )

ADD_LIBRARY ( pbrt STATIC
  ${PBRT_YACC_LEX_SOURCE}
  ${PBRT_CORE_SOURCE}
  ${PBRT_CORE_HEADERS}
  ${PBRT_SOURCE}
  )
ADD_SANITIZERS ( pbrt )

IF (WIN32)
  # Avoid a name clash when building on Visual Studio
  SET_TARGET_PROPERTIES ( pbrt
    PROPERTIES
    OUTPUT_NAME libpbrt
    )
ENDIF()

# Main renderer
ADD_EXECUTABLE ( pbrt_exe
  src/main/pbrt.cpp
  )
ADD_SANITIZERS ( pbrt_exe )

SET_TARGET_PROPERTIES ( pbrt_exe
  PROPERTIES
  OUTPUT_NAME pbrt
  )

# Create a custom set of library list based on the build type
option(ENABLE_OPENVDB_SUPPORT "Build with openvdb support" ON)

# also set some compile definitions for libpbrt and pbrt.exec_program
# the following flag is set so that if openvdb was compiled with that flag, this is required
# to correctly link them with pbrt
option(OPENVDB_ABI3 "Use openvdb that is abi 3 compatible" OFF)

if (OPENVDB_ABI3)
	target_compile_definitions(pbrt PRIVATE OPENVDB_3_ABI_COMPATIBLE)
	target_compile_definitions(pbrt_exe PRIVATE OPENVDB_3_ABI_COMPATIBLE)
endif(OPENVDB_ABI3)

# the following flag is set in the code
# we appropriately set it through cmake
target_compile_definitions(pbrt PRIVATE PBRT_FLOAT_AS_DOUBLE)
target_compile_definitions(pbrt PRIVATE ENABLE_OPENVDB_SUPPORT)
#else()
#set (REQUIRED_LIBS 
#  pbrt
#  ${CMAKE_THREAD_LIBS_INIT}
#  ${OPENEXR_LIBS}
#  glog
#)
##endif()

option(ENABLE_FIELD3D_SUPPORT "Build with Field3d support" ON)
target_compile_definitions(pbrt PRIVATE ENABLE_FIELD3D_SUPPORT)
target_compile_definitions(pbrt PRIVATE FIELD3D_STATIC)
target_compile_definitions(pbrt_exe PRIVATE FIELD3D_STATIC)
target_compile_definitions(pbrt_exe PRIVATE ENABLE_FIELD3D_SUPPORT)

set(FIELD3D_LIBRARY 
	optimized ${FIELD3D_LIBRARY_RELEASE} 
	debug ${FIELD3D_LIBRARY_DEBUG})

set(HDF5_LIBRARIES 
	optimized ${HDF5_LIBRARY_RELEASE} ${HDF5_CPP_LIBRARY_RELEASE} ${SZIP_LIBRARY_RELEASE} ${ZLIB_LIBRARY_RELEASE} 
	debug ${HDF5_LIBRARY_DEBUG} ${HDF5_CPP_LIBRARY_DEBUG} ${SZIP_LIBRARY_DEBUG} ${ZLIB_LIBRARY_DEBUG})

if(WIN32)
set (REQUIRED_LIBS
  pbrt
  ${CMAKE_THREAD_LIBS_INIT}
  ${OPENEXR_LIBS}
  glog
  ${FIELD3D_LIBRARY}
  ${Boost_LIBRARIES}
  ${HDF5_LIBRARIES}
)
else()
# we need two additional libraries dl and rt for successful linking
set (REQUIRED_LIBS 
  pbrt
  ${CMAKE_THREAD_LIBS_INIT}
  dl
  glog
  ${Boost_LIBRARIES}
  ${FIELD3D_LIBRARY}
  ${HDF5_LIBRARIES}
  ${Ilmbase_HALF_LIBRARY}
  ${OPENEXR_LIBS}
  ${OPENVDB_LIBRARY}
  ${Tbb_TBB_LIBRARY}
  ${Tbb_TBBMALLOC_LIBRARY}
)
endif() # end of win32 conditional check
message("Important Libs--->")
message("\t-Field3D:" ${FIELD3D_LIBRARY})
message("\t-OPENEXR:" ${OPENEXR_LIBS})
message("\t-OpenVDB:" ${OPENVDB_LIBRARY})
message("\t\t include dir:" ${OPENVDB_INCLUDE_DIR})
message("\t-TBB:" ${Tbb_TBB_LIBRARY})
message("\t-TBB Malloc:" ${Tbb_TBBMALLOC_LIBRARY})

TARGET_LINK_LIBRARIES (pbrt_exe ${REQUIRED_LIBS})

# Tools
#ADD_EXECUTABLE ( bsdftest
#  src/tools/bsdftest.cpp
#  )
#ADD_SANITIZERS ( bsdftest )

#ADD_EXECUTABLE ( imgtool
#  src/tools/imgtool.cpp
#  )
#ADD_SANITIZERS ( imgtool )

#ADD_EXECUTABLE ( obj2pbrt
#  src/tools/obj2pbrt.cpp
#  )
#ADD_SANITIZERS ( obj2pbrt )

#ADD_EXECUTABLE ( cyhair2pbrt
#  src/tools/cyhair2pbrt.cpp
#  )
#ADD_SANITIZERS ( cyhair2pbrt )

#TARGET_LINK_LIBRARIES ( bsdftest
#  ${REQUIRED_LIBS}
#  )

#TARGET_LINK_LIBRARIES ( imgtool
#  ${REQUIRED_LIBS}
#  )

# Unit test

#FILE ( GLOB PBRT_TEST_SOURCE
#  src/tests/*.cpp
#  src/tests/gtest/*.cc
#  )

#ADD_EXECUTABLE ( pbrt_test
#  ${PBRT_TEST_SOURCE}
#  )
#ADD_SANITIZERS ( pbrt_test )

#TARGET_LINK_LIBRARIES ( pbrt_test
#  ${REQUIRED_LIBS}
#)

# Installation

#INSTALL ( TARGETS
#  pbrt_exe
  #bsdftest
  #imgtool
  #  obj2pbrt
  #cyhair2pbrt
  #  DESTINATION
  #  bin
  #  )

INSTALL ( TARGETS
  pbrt
  DESTINATION
  lib
)
