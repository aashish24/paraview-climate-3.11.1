PROJECT(vtkVisItAVTAlgorithms)

SET(KIT VisItAVTAlgorithms)
SET(KIT_UPCASE "VISITAVTALGORITHMS")

SET(Kit_SRCS
vtkAvtFileFormatAlgorithm.h
vtkAvtFileFormatAlgorithm.cxx
vtkAvtMTSDFileFormatAlgorithm.h
vtkAvtMTSDFileFormatAlgorithm.cxx
vtkAvtMTMDFileFormatAlgorithm.h
vtkAvtMTMDFileFormatAlgorithm.cxx
vtkAvtSTSDFileFormatAlgorithm.h
vtkAvtSTSDFileFormatAlgorithm.cxx
vtkAvtSTMDFileFormatAlgorithm.h
vtkAvtSTMDFileFormatAlgorithm.cxx
)

#setup the kit so we get client server bindings  
SET(VTK_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")
INCLUDE(${VTK_CMAKE_DIR}/vtkExportKit.cmake)
VTK_EXPORT_KIT("${KIT}" "${KIT_UPCASE}" "${Kit_SRCS}")

#########################################################################  
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/vtkAvtAlgorithmsExport.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/vtkAvtAlgorithmsExport.h")
 set(dynamicHeaders "${dynamicHeaders};${CMAKE_CURRENT_BINARY_DIR}/vtkAvtAlgorithmsExport.h")

INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}
${CMAKE_CURRENT_BINARY_DIR}
${AVT_DATABASE_SOURCE_DIR}/Database
${AVT_DATABASE_SOURCE_DIR}/Ghost
${AVT_DATABASE_SOURCE_DIR}/Formats
${VISIT_COMMON_INCLUDES}
${AVT_DBATTS_SOURCE_DIR}/MetaData
${AVT_DBATTS_SOURCE_DIR}/SIL
${VISIT_SOURCE_DIR}/avt/Math
${VISIT_SOURCE_DIR}/avt/VisWindow/VisWindow
${AVT_PIPELINE_SOURCE_DIR}/AbstractFilters
${AVT_PIPELINE_SOURCE_DIR}/Data
${AVT_PIPELINE_SOURCE_DIR}/Pipeline
${AVT_PIPELINE_SOURCE_DIR}/Sinks
${AVT_PIPELINE_SOURCE_DIR}/Sources
${VISIT_SOURCE_DIR}/visit_vtk/full
${VISIT_SOURCE_DIR}/visit_vtk/lightweight
${VTK_INCLUDE_DIRS}
)
VTK_ADD_LIBRARY(vtkVisItAVTAlgorithms ${Kit_SRCS})
IF(VISIT_PARALLEL)    
  TARGET_LINK_LIBRARIES(vtkVisItAVTAlgorithms visitcommon avtdatabase_par avtpipeline_par visit_vtk avtivp)    
ELSE(VISIT_PARALLEL)	
  TARGET_LINK_LIBRARIES(vtkVisItAVTAlgorithms visitcommon avtdatabase_ser avtpipeline_ser visit_vtk avtivp)
ENDIF(VISIT_PARALLEL)

#Wrap for Python
#-------------------------------------------------------------------------------------------------
IF (PARAVIEW_ENABLE_PYTHON)  
  SET(Kit_PYTHON_EXTRA_SRCS)
  SET(KIT_PYTHON_LIBS vtkVisItAVTAlgorithmsPythonD vtkIOPythonD vtkVisItAVTAlgorithmsCS)
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
VISIT_INSTALL_TARGETS(vtkVisItAVTAlgorithms)