MACRO(VISIT_INCLUDE_READER_SOURCE name state)
option(VISIT_BUILD_READER_${name} "Build VisIt Reader ${name}" ${state})
if (VISIT_BUILD_READER_${name})
  ADD_SUBDIRECTORY(${name})
endif()
ENDMACRO(VISIT_INCLUDE_READER_SOURCE name)

PROJECT(vtkVisItDatabases)

#these are used by visit_include_reader
set(VISIT_DATABASE_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(VISIT_SERVER_SOURCES CACHE INTERNAL "vtk classes to wrap for client server" FORCE)
set(VISIT_DB_SOURCES "" CACHE INTERNAL "visit sources for readers" FORCE) 
set(VISIT_DB_INC_DIRS CACHE INTERNAL "include directories" FORCE)

set(DEFAULT_BRIDGE_READERS
  ANALYZE
  ANSYS  
  AUXFile
  BOV
  CEAucd
  Chombo  
  CMAT
  CTRL
  Curve2D
  DDCMD
  Dyna3D
  Enzo
  ExtrudedVol
  FLASH  
  Fluent
  Gadget
  GGCM
  GTC
  GULP
  H5Nimrod
  Image
  LAMMPS
  Lines
  M3D
  M3DC1
  Miranda
  MM5
  NASTRAN
  NETCDF
  Nek5000
  OpenFOAM
  OVERFLOW
  paraDIS
  PATRAN
  PFLOTRAN
  Pixie
  PlainText
  PLOT2D
  PLOT3D
  Point3D 
  ProteinDataBank
  RAW
  SAMRAI
  SAR
  SAS  
  Spheral
  Tecplot  
  Tetrad
  TFT
  TSurf
  UNIC
  VASP
  Velodyne
  Vs
  Xmdv
  XYZ
  )
if(NOT WIN32)
  list(APPEND DEFAULT_BRIDGE_READERS
    Claw
    )
endif(NOT WIN32)    
    
#add all the default readers to the library    
foreach(db ${DEFAULT_BRIDGE_READERS})
  add_subdirectory(${db})
endforeach(db)

ADD_DEFINITIONS("-DRC_CPP_VISIT_BUILD")
  
VISIT_INCLUDE_READER_SOURCE(Silo OFF)
if(VISIT_BUILD_READER_Silo)
  if(NOT HAVE_LIBSILO)
    #if we currently don't have silo setup make it required
    find_package(SILO REQUIRED)
  endif(NOT HAVE_LIBSILO)
  if(NOT SILO_FOUND)
    message(FATAL_ERROR "Unable to locate Silo library.")
  endif(NOT SILO_FOUND)
  set(VISITBRIDGE_USE_SILO ${SILO_FOUND} PARENT_SCOPE)
endif(VISIT_BUILD_READER_Silo)


VISIT_INCLUDE_READER_SOURCE(CGNS OFF)
if(VISIT_BUILD_READER_CGNS)
  if(NOT HAVE_LIBCGNS)
    #if we currently don't have silo setup make it required
    find_package(CGNS REQUIRED)
  endif(NOT HAVE_LIBCGNS)
  if(NOT CGNS_FOUND)
    message(FATAL_ERROR "Unable to locate CGNS library.")
  endif(NOT CGNS_FOUND)
  set(VISITBRIDGE_USE_CGNS ${CGNS_FOUND} PARENT_SCOPE)
endif(VISIT_BUILD_READER_CGNS)

VISIT_INCLUDE_READER_SOURCE(Mili OFF)
if(VISIT_BUILD_READER_Mili)
  if(NOT HAVE_LIBMILI)
    #if we currently don't have silo setup make it required
    find_package(MILI REQUIRED)
  endif(NOT HAVE_LIBMILI)
  if(NOT MILI_FOUND)
    message(FATAL_ERROR "Unable to locate Mili library.")
  endif(NOT MILI_FOUND)
  set(VISITBRIDGE_USE_MILI ${MILI_FOUND} PARENT_SCOPE)  
endif(VISIT_BUILD_READER_Mili)

SET(KIT VisItDatabases)
SET(KIT_UPCASE "VISITDATABASES")

#include all the db directores before
#the visit files
include_directories(BEFORE
  ${VISIT_DB_INC_DIRS}
  )

include_directories(
  ${VISITBRIDGE_INCLUDE_DIRS}
  )

if(HAVE_LIBSILO AND VISIT_BUILD_READER_Silo)
  include_directories(${SILO_INCLUDE_DIR})
endif(HAVE_LIBSILO AND VISIT_BUILD_READER_Silo)

if(HAVE_LIBCGNS AND VISIT_BUILD_READER_CGNS)
  include_directories(${CGNS_INCLUDE_DIR})
endif(HAVE_LIBCGNS AND VISIT_BUILD_READER_CGNS)

if(HAVE_LIBMILI AND VISIT_BUILD_READER_Mili)
  include_directories(${MILI_INCLUDE_DIR})
endif(HAVE_LIBMILI AND VISIT_BUILD_READER_Mili)
  
SET(Kit_SRCS
  ${VISIT_SERVER_SOURCES}
  )

SET_SOURCE_FILES_PROPERTIES(
  ${VISIT_DB_SOURCES}
  WRAP_EXCLUDE
  )

#setup the kit so we get client server bindings  
SET(VTK_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")
INCLUDE(${VTK_CMAKE_DIR}/vtkExportKit.cmake)
VTK_EXPORT_KIT("${KIT}" "${KIT_UPCASE}" "${Kit_SRCS}")
  
#build all the readers into a single library
vtk_add_library(vtk${KIT} ${Kit_SRCS} ${VISIT_DB_SOURCES})

target_link_libraries(vtk${KIT} vtkVisItAVTAlgorithms vtkNetCDF) 
if(PARAVIEW_HDF5_LIBRARIES)
  target_link_libraries(vtk${KIT} ${PARAVIEW_HDF5_LIBRARIES})
endif(PARAVIEW_HDF5_LIBRARIES)  

if(HAVE_LIBSILO AND VISIT_BUILD_READER_Silo)
  target_link_libraries(vtk${KIT} ${SILO_LIBRARIES})
endif(HAVE_LIBSILO AND VISIT_BUILD_READER_Silo)

if(HAVE_LIBCGNS AND VISIT_BUILD_READER_CGNS)  
  target_link_libraries(vtk${KIT} ${CGNS_LIBRARIES})
endif(HAVE_LIBCGNS AND VISIT_BUILD_READER_CGNS)

if(HAVE_LIBMILI AND VISIT_BUILD_READER_Mili)  
  target_link_libraries(vtk${KIT} ${MILI_LIBRARIES})
endif(HAVE_LIBMILI AND VISIT_BUILD_READER_Mili)

#Wrap for Python
#-------------------------------------------------------------------------------------------------
IF (PARAVIEW_ENABLE_PYTHON)  
  SET(Kit_PYTHON_EXTRA_SRCS)
  SET(KIT_PYTHON_LIBS vtk${KIT}PythonD vtk${KIT}CS
    vtkVisItAVTAlgorithmsPythonD vtkVisItAVTAlgorithmsCS)
  # Tell vtkWrapPython.cmake to set VTK_PYTHON_LIBRARIES for us.
  SET(VTK_WRAP_PYTHON_FIND_LIBS 1)
  SET(VTK_WRAP_INCLUDE_DIRS ${VTK_INCLUDE_DIR})
  INCLUDE("${VTK_CMAKE_DIR}/vtkWrapPython.cmake")
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  SET(KIT_PYTHON_DEPS)
  SET(VTK_INSTALL_NO_LIBRARIES 1)
  IF(VTKPythonWrapping_INSTALL_BIN_DIR)
    SET(VTK_INSTALL_NO_LIBRARIES)
  ENDIF(VTKPythonWrapping_INSTALL_BIN_DIR)

  SET(VTK_INSTALL_LIB_DIR "${VTKPythonWrapping_INSTALL_LIB_DIR}")
  SET(VTK_INSTALL_BIN_DIR "${VTKPythonWrapping_INSTALL_BIN_DIR}")
  SET(VTK_INSTALL_LIB_DIR_CM24 "${VTKPythonWrapping_INSTALL_LIB_DIR}")
  SET(VTK_INSTALL_BIN_DIR_CM24 "${VTKPythonWrapping_INSTALL_BIN_DIR}")

  # include python wrapping from VTK
  INCLUDE(KitCommonPythonWrapBlock)
ENDIF (PARAVIEW_ENABLE_PYTHON)

#Install
#-------------------------------------------------------------------------------------------------
VISIT_INSTALL_TARGETS(vtk${KIT})
