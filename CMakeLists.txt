# ==========================
# BamTools CMakeLists.txt
# (c) 2010 Derek Barnett
#
# top-level
# ==========================

# set project name
project( BamTools )

# Cmake requirements
cmake_minimum_required( VERSION 2.6.4 )

# Force the build directory to be different from source directory
macro( ENSURE_OUT_OF_SOURCE_BUILD MSG )
    string( COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource )
    get_filename_component( PARENTDIR ${CMAKE_SOURCE_DIR} PATH )
    string( COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${PARENTDIR}" insourcesubdir )
    IF( insource OR insourcesubdir )
        message( FATAL_ERROR "${MSG}" )
    ENDIF( insource OR insourcesubdir )
endmacro( ENSURE_OUT_OF_SOURCE_BUILD )

ensure_out_of_source_build( "
  ${PROJECT_NAME} requires an out of source build.
  $ mkdir build 
  $ cd build
  $ cmake ..
  $ make 
(or the Windows equivalent)\n" )

# set BamTools version information
set( BamTools_VERSION_MAJOR 2 )
set( BamTools_VERSION_MINOR 5 )
set( BamTools_VERSION_BUILD 6 )

# set our library and executable destination dirs
set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )
set( LIBRARY_OUTPUT_PATH    "${CMAKE_SOURCE_DIR}/lib" )

#set(CMAKE_BUILD_TYPE "Debug")
# to use the release mode
# cmake . -DCMAKE_BUILD_TYPE=Release

# define compiler flags for all code
set( CMAKE_BUILD_TYPE Release )
# old setting std=c++98
set( CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -Wall ${CMAKE_CXX_FLAGS_RELEASE}" )
set( CMAKE_CXX_FLAGS_DEBUG "-g -std=c++11 -Wall ${CMAKE_CXX_FLAGS_DEBUG}" )
# add additional compiler flags use add_definitions()
add_definitions( -D_FILE_OFFSET_BITS=64 )
message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")

find_package(ZLIB)
if (ZLIB_FOUND)
	include_directories(${ZLIB_INDLUDE_DIRS})
endif (ZLIB_FOUND)

# -----------------------------------------------
# handle platform-/environment-specific defines

# If planning to run in Node.js environment, run:
# cmake -DEnableNodeJS=true
if( EnableNodeJS )
    add_definitions( -DSYSTEM_NODEJS=1 )
endif()

# If running on SunOS
if( "${CMAKE_SYSTEM_NAME}" MATCHES "SunOS" )
    add_definitions( -DSUN_OS )
endif()

# -------------------------------------------

# add our includes root path
include_directories( src )

# list subdirectories to build in
add_subdirectory( src )
