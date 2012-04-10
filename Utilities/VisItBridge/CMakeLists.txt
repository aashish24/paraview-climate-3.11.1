cmake_minimum_required(VERSION 2.6.3)

project(VisItBridgePlugin)

#Setup the version since some files need this
SET(VISIT_VERSION "2.0.0")
    
IF (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
ENDIF (COMMAND cmake_policy)

#setup the INSTALL paths for the kits
set(VTK_INSTALL_PACKAGE_DIR_CM24 ${PV_INSTALL_LIB_DIR})  
set(VTK_INSTALL_PACKAGE_DIR ${PV_INSTALL_LIB_DIR})

#disable compiler warnings from the bridge
option(VISIT_DISABLE_COMPILER_WARNINGS "Disable compiler warnings" ON)
mark_as_advanced(VISIT_DISABLE_COMPILER_WARNINGS)
if(VISIT_DISABLE_COMPILER_WARNINGS)
  if(WIN32)
    if (MSVC)    
      string(REGEX REPLACE "(^| )([/-])W[0-9]( |$)" " "
          CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
      set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} /w")
      
      string(REGEX REPLACE "(^| )([/-])W[0-9]( |$)" " "
          CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /w")
    endif(MSVC)  
  else(WIN32)  
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    set (CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -w")
  endif(WIN32)
endif(VISIT_DISABLE_COMPILER_WARNINGS)

#set up some vars we need to compile
set(VISIT_STATIC)
if (BUILD_SHARED_LIBS)
  set(VISIT_STATIC 0)  
  add_definitions(-DVISIT_BUILD_SHARED_LIBS)
else(VISIT_STATIC)
  set(VISIT_STATIC 1)
  add_definitions(-DVISIT_STATIC -DGLEW_STATIC)
endif()

set(VISIT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(VISIT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(VISIT_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMake )

#include the visit cmake directory on the cmake search path
SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${VISIT_CMAKE_DIR})

#include the visit install target and plugin function
include(${VISIT_CMAKE_DIR}/VisItCommon.cmake)

#need toi grab some helper macros from paraview plugin 
INCLUDE(ParaViewPlugins)
 
#set up MPI
set(VISIT_PARALLEL ${VTK_USE_MPI})
if(VTK_USE_MPI)
  INCLUDE(FindMPI)
  INCLUDE_DIRECTORIES(  
    ${MPI_INCLUDE_PATH}
    )  
  set(VISIT_PARALLEL_LIBS ${MPI_LIBRARY})
  if(MPI_EXTRA_LIBRARY)
  set(VISIT_PARALLEL_LIBS ${VISIT_PARALLEL_LIBS} ${MPI_EXTRA_LIBRARY})
  endif(MPI_EXTRA_LIBRARY)
endif(VTK_USE_MPI)


#block out most of the warnings in visit on windows
IF (WIN32)
  ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
  IF (MSVC_VERSION GREATER 1400)
      ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
      ADD_DEFINITIONS(-D_SCL_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE)
      ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
  ENDIF (MSVC_VERSION GREATER 1400)
ENDIF(WIN32)  

#-----------------------------------------------------------------------------
# Setup ParaView and Common includes before the visit-config.h so that
# we can use utilities like HDF5 from ParaView
#-----------------------------------------------------------------------------


# Set up easy to use includes for the common directory
set(VISIT_COMMON_INCLUDES
    ${CMAKE_CURRENT_BINARY_DIR}  
    ${CMAKE_CURRENT_SOURCE_DIR}/common/Exceptions/Database
    ${CMAKE_CURRENT_SOURCE_DIR}/common/Exceptions/Pipeline
    ${CMAKE_CURRENT_SOURCE_DIR}/common/Exceptions/Plotter
    ${CMAKE_CURRENT_SOURCE_DIR}/common/comm
    ${CMAKE_CURRENT_SOURCE_DIR}/common/expr
    ${CMAKE_CURRENT_SOURCE_DIR}/common/icons
    ${CMAKE_CURRENT_SOURCE_DIR}/common/misc  
    ${CMAKE_CURRENT_SOURCE_DIR}/common/parser
    ${CMAKE_CURRENT_SOURCE_DIR}/common/plugin
    ${CMAKE_CURRENT_SOURCE_DIR}/common/proxybase
    ${CMAKE_CURRENT_SOURCE_DIR}/common/state
    ${CMAKE_CURRENT_SOURCE_DIR}/common/utility
    ${CMAKE_CURRENT_SOURCE_DIR}/common/common/misc
    ${CMAKE_CURRENT_SOURCE_DIR}/common/common/plugin
    ${CMAKE_CURRENT_SOURCE_DIR}/common/common/state
    ${CMAKE_CURRENT_SOURCE_DIR}/common/common/utility
  )
  
#watch out, you need to make sure common/parser is always in front of
# python2.X includes
INCLUDE_DIRECTORIES(BEFORE ${VISIT_COMMON_INCLUDES})
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})

SET(VTK_BINARY_DIR ${PARAVIEW_VTK_DIR} )   

#-----------------------------------------------------------------------------
# Setup Vars for visit-config.h
#-----------------------------------------------------------------------------


set(VISIT_DBIO_ONLY ON) #Build only visitconvert and engine plugins
if(VISIT_DBIO_ONLY)
    add_definitions(-DDBIO_ONLY)    
endif(VISIT_DBIO_ONLY)   

#Check to see if ofstreams rdbuf is public. If it is NOT public set NO_SETBUF
TRY_COMPILE(VISIT_COMPILER_FSTREAM_WORKAROUND
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp
        ${VISIT_CMAKE_DIR}/testFStream.cxx)
IF (NOT VISIT_COMPILER_FSTREAM_WORKAROUND)
   set(NO_SETBUF 1)
ENDIF (NOT VISIT_COMPILER_FSTREAM_WORKAROUND)  

# Set the slash characters based on the platform
IF(WIN32)
    SET(VISIT_SLASH_CHAR   "'\\\\'")
    SET(VISIT_SLASH_STRING "\"\\\\\"")
ELSE(WIN32)
    SET(VISIT_SLASH_CHAR   "'/'")
    SET(VISIT_SLASH_STRING "\"/\"")
ENDIF(WIN32)

# Check for plugin extension
IF(VISIT_STATIC)
  IF(WIN32)
    SET(VISIT_PLUGIN_EXTENSION   ".dll")
  ELSE(WIN32)
    SET(VISIT_PLUGIN_EXTENSION   ".lib")    
  ENDIF(WIN32)
ELSE(VISIT_STATIC)
    IF(WIN32)
        SET(VISIT_PLUGIN_EXTENSION   ".dll")
    ELSE(WIN32)
        IF(APPLE)
            SET(VISIT_PLUGIN_EXTENSION   ".dylib")
        ELSE(APPLE)
            SET(VISIT_PLUGIN_EXTENSION   ".so")
        ENDIF(APPLE)
    ENDIF(WIN32)
ENDIF(VISIT_STATIC)

#-----------------------------------------------------------------------------
# Setup lib settings
#-----------------------------------------------------------------------------

# Set up boost (interval) library
find_package(Boost 1.40 REQUIRED) #only need header libraries
if (Boost_Found)  
  set(HAVE_BILIB 1)
else (Boost_Found)
  set(HAVE_BILIB 0)
endif(Boost_Found)


#setup non third party vtk utilities
set(HAVE_LIBHDF5 ${PARAVIEW_HDF5_LIBRARIES})
set(HAVE_LIBNETCDF 1)
set(HAVE_NETCDF_H 1)
set(HAVE_LIBEXODUSII 1)

# Setup SILO
find_package(SILO QUIET)
if(SILO_FOUND)
  set(HAVE_LIBSILO ${SILO_FOUND})
  set(VISITBRIDGE_USE_SILO ${SILO_FOUND} PARENT_SCOPE)
endif(SILO_FOUND)

# Setup CGNS
find_package(CGNS QUIET)
if(CGNS_FOUND)
  set(HAVE_LIBCGNS ${CGNS_FOUND})
  set(VISITBRIDGE_USE_CGNS ${CGNS_FOUND} PARENT_SCOPE)
endif(CGNS_FOUND)

# Setup Mili
find_package(MILI QUIET)
if(MILI_FOUND)
  set(HAVE_LIBMILI ${MILI_FOUND})
  set(VISITBRIDGE_USE_MILI ${MILI_FOUND} PARENT_SCOPE)
endif(MILI_FOUND)

#setup zlib
include_directories(${VTK_ZLIB_INCLUDE_DIRS})

#setup netcdf
VISIT_VTK_THIRD_PARTY_INCLUDE(NETCDF vtknetcdf)

#setup png
VISIT_VTK_THIRD_PARTY_INCLUDE(PNG vtkpng)


#-----------------------------------------------------------------------------
# Detect packages here. We could probably write macros that we can include from
# elsewhere for this. 
#-----------------------------------------------------------------------------
INCLUDE(CheckIncludeFiles)
INCLUDE(CMakeBackwardCompatibilityC)
INCLUDE(CMakeBackwardCompatibilityCXX)
INCLUDE(CheckTypeSize)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckSymbolExists)
INCLUDE(TestBigEndian)
INCLUDE(FindOpenGL)

CHECK_INCLUDE_FILES (fcntl.h     HAVE_FCNTL_H)
CHECK_INCLUDE_FILES (inttypes.h  HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES (malloc.h    HAVE_MALLOC_H)
CHECK_INCLUDE_FILES (limits.h    HAVE_LIMITS_H)
CHECK_INCLUDE_FILES (memory.h    HAVE_MEMORY_H)
CHECK_INCLUDE_FILES (stdint.h    HAVE_STDINT_H)
CHECK_INCLUDE_FILES (stdlib.h    HAVE_STDLIB_H)
CHECK_INCLUDE_FILES (strings.h   HAVE_STRINGS_H)
CHECK_INCLUDE_FILES (string.h    HAVE_STRING_H)
CHECK_INCLUDE_FILES (sys/time.h  HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES (sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES (sys/stat.h  HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES (unistd.h    HAVE_UNISTD_H)
CHECK_INCLUDE_FILES (stdbool.h   HAVE_STDBOOL_H)

# Check for type sizes, endian
SET(SIZEOF_BOOLEAN              ${CMAKE_SIZEOF_BOOLEAN})
SET(SIZEOF_CHAR                 ${CMAKE_SIZEOF_CHAR})
SET(SIZEOF_DOUBLE               ${CMAKE_SIZEOF_DOUBLE})
SET(SIZEOF_FLOAT                ${CMAKE_SIZEOF_FLOAT})
SET(SIZEOF_INT                  ${CMAKE_SIZEOF_INT})
SET(SIZEOF_LONG                 ${CMAKE_SIZEOF_LONG})
SET(SIZEOF_LONG_DOUBLE          ${CMAKE_SIZEOF_LONG_DOUBLE})
SET(SIZEOF_LONG_FLOAT           ${CMAKE_SIZEOF_LONG_FLOAT})
SET(SIZEOF_LONG_LONG            ${CMAKE_SIZEOF_LONG_LONG})
SET(SIZEOF_SHORT                ${CMAKE_SIZEOF_SHORT})
SET(SIZEOF_UNSIGNED_CHAR        ${CMAKE_SIZEOF_UNSIGNED_CHAR})
SET(SIZEOF_UNSIGNED_INT         ${CMAKE_SIZEOF_UNSIGNED_INT})
SET(SIZEOF_UNSIGNED_LONG        ${CMAKE_SIZEOF_UNSIGNED_LONG})
SET(SIZEOF_UNSIGNED_LONG_LONG   ${CMAKE_SIZEOF_UNSIGNED_LONG_LONG})
SET(SIZEOF_UNSIGNED_SHORT       ${CMAKE_SIZEOF_UNSIGNED_SHORT})
SET(SIZEOF_VOID_P               ${CMAKE_SIZEOF_VOID_P})
SET(CMAKE_REQUIRED_DEFINITIONS -D_LARGEFILE64_SOURCE)
CHECK_TYPE_SIZE("off64_t" SIZEOF_OFF64_T)
TEST_BIG_ENDIAN(WORDS_BIGENDIAN)


# manually check for socklen_t as CHECK_SYMBOL_EXISTS
# doesn't appear to work on linux (at least) 
TRY_COMPILE(HAVE_SOCKLEN_T
    ${CMAKE_CURRENT_BINARY_DIR}
    ${VISIT_SOURCE_DIR}/CMake/TestSocklenT.c
    OUTPUT_VARIABLE SLT
)
IF (HAVE_SOCKLEN_T)
    SET(HAVE_SOCKLEN_T 1 CACHE INTERNAL "support for socklen_t")
ELSE(HAVE_SOCKLEN_T)
    SET(HAVE_SOCKLEN_T 0 CACHE INTERNAL "support for socklen_t")
ENDIF (HAVE_SOCKLEN_T)

#-----------------------------------------------------------------------------
# Create visit-config.h
#-----------------------------------------------------------------------------

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/visit-cmake.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/visit-config.h @ONLY IMMEDIATE)  

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/common/misc/ConfigureInfo.C.in
               ${CMAKE_CURRENT_BINARY_DIR}/common/misc/ConfigureInfo.C @ONLY IMMEDIATE)

               
add_subdirectory(common)
add_subdirectory(visit_vtk)
add_subdirectory(avt)
add_subdirectory(AvtAlgorithms)

set(VISITBRIDGE_INCLUDE_DIRS    
  ${vtkVisItAVTAlgorithms_BINARY_DIR}
  ${vtkVisItAVTAlgorithms_SOURCE_DIR}
  ${VISIT_COMMON_INCLUDES}
  ${VISIT_SOURCE_DIR}/include
  ${AVT_DATABASE_SOURCE_DIR}/Database
  ${AVT_DATABASE_SOURCE_DIR}/Ghost
  ${AVT_DATABASE_SOURCE_DIR}/Formats
  ${AVT_DBATTS_SOURCE_DIR}/MetaData
  ${AVT_DBATTS_SOURCE_DIR}/SIL
  ${VISIT_SOURCE_DIR}/avt/Math
  ${VISIT_SOURCE_DIR}/avt/VisWindow/VisWindow
  ${VISIT_SOURCE_DIR}/avt/IVP
  ${AVT_PIPELINE_SOURCE_DIR}/AbstractFilters
  ${AVT_PIPELINE_SOURCE_DIR}/Data
  ${AVT_PIPELINE_SOURCE_DIR}/Pipeline
  ${AVT_PIPELINE_SOURCE_DIR}/Sinks
  ${AVT_PIPELINE_SOURCE_DIR}/Sources
  ${VISIT_SOURCE_DIR}/visit_vtk/full
  ${VISIT_SOURCE_DIR}/visit_vtk/lightweight
  ${Boost_INCLUDE_DIR}
  ${VTK_INCLUDE_DIR}
  )
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/VisItBridgeUse.cmake.in
               ${CMAKE_CURRENT_BINARY_DIR}/VisItBridgeUse.cmake @ONLY IMMEDIATE)



add_subdirectory(databases)


#######################################################################
# Update if we found SILO,CGNS,MILI
  set(VISITBRIDGE_USE_MILI ${VISITBRIDGE_USE_MILI} PARENT_SCOPE)  
  set(VISITBRIDGE_USE_CGNS ${VISITBRIDGE_USE_CGNS} PARENT_SCOPE)
  set(VISITBRIDGE_USE_SILO ${VISITBRIDGE_USE_SILO} PARENT_SCOPE)

#######################################################################
# INSTALL
IF (NOT PV_INSTALL_NO_DEVELOPMENT)
  SET(BUILD_TREE ${ParaView_BINARY_DIR})
  
  #install
  SET(cmakeDevFiles
    ${VISIT_CMAKE_DIR}/VisItCommon.cmake
    ${VISIT_CMAKE_DIR}/VisItExport.h.in
    ${VISIT_CMAKE_DIR}/VisItGUI.xml.in
    ${VISIT_CMAKE_DIR}/VisItMT.h.in
    ${VISIT_CMAKE_DIR}/VisItMT.cxx.in
    ${VISIT_CMAKE_DIR}/VisItSM.xml.in
    ${VISIT_CMAKE_DIR}/VisItST.h.in
    ${VISIT_CMAKE_DIR}/VisItST.cxx.in)
      
  INSTALL(
    FILES ${cmakeDevFiles} 
    DESTINATION ${PV_INSTALL_CMAKE_DIR}
    COMPONENT Development)
      
  SET(VISITBRIDGE_INCLUDE_DIRS "\${ParaView_INSTALL_PREFIX}/${PV_INSTALL_INCLUDE_DIR}")
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/VisItBridgeUse.cmake.in
                "${BUILD_TREE}/CMake/tmp/VisItBridgeUse.cmake" @ONLY IMMEDIATE)  
  INSTALL(
        FILES "${BUILD_TREE}/CMake/tmp/VisItBridgeUse.cmake"          
        DESTINATION ${PV_INSTALL_CMAKE_DIR}
        COMPONENT Development)  
ENDIF (NOT PV_INSTALL_NO_DEVELOPMENT)
