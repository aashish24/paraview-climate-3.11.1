# PVVTKExtensions should contain classes that extend VTK. This can only depend
# on VTK and PVCommon.

PROJECT(PVVTKExtensions)

INCLUDE_REGULAR_EXPRESSION("^vtk.*$")

# needed by vtkExportKit.cmake
SET(VTK_INSTALL_PACKAGE_DIR_CM24 ${PV_INSTALL_LIB_DIR})
SET(VTK_INSTALL_PACKAGE_DIR ${PV_INSTALL_LIB_DIR})

SET(KIT PVVTKExtensions)
SET(KIT_UPCASE "PVVTKEXTENSIONS")
IF(NOT WIN32)
  SET(CMAKE_C_FLAGS "${CMAKE_ANSI_CFLAGS} ${CMAKE_C_FLAGS}")
  ADD_DEFINITIONS(-D_HPUX_SOURCE)
ENDIF(NOT WIN32)

SET (Kit_SRCS
  vtkAbstractPOPReader.cxx
  vtkAMRDualClip.cxx
  vtkAMRDualContour.cxx
  vtkAMRDualGridHelper.cxx
  vtkAnimationPlayer.cxx
  vtkAppendArcLength.cxx
  vtkAppendRectilinearGrid.cxx
  vtkAttributeDataReductionFilter.cxx
  vtkAttributeDataToTableFilter.cxx
  vtkBlockDeliveryPreprocessor.cxx
  vtkBSPCutsGenerator.cxx
  vtkCacheSizeKeeper.cxx
  vtkCameraInterpolator2.cxx
  vtkCameraManipulator.cxx
  vtkCameraManipulatorGUIHelper.cxx
  vtkCellIntegrator.cxx
  vtkCleanArrays.cxx
  vtkCleanUnstructuredGrid.cxx
  vtkCleanUnstructuredGridCells.cxx
  vtkCompositeAnimationPlayer.cxx
  vtkCompositeDataToUnstructuredGridFilter.cxx
  vtkContextNamedOptions.cxx
  vtkCSVExporter.cxx
  vtkCSVWriter.cxx
  vtkDataSetToRectilinearGrid.cxx
  vtkEnzoReader.cxx
  vtkEquivalenceSet.cxx
  vtkExodusFileSeriesReader.cxx
  vtkExtractHistogram.cxx
  vtkExtractScatterPlot.cxx
  vtkFileSequenceParser.cxx
  vtkFileSeriesReader.cxx
  vtkFileSeriesWriter.cxx
  vtkFlashContour.cxx
  vtkFlashReader.cxx
  vtkGridConnectivity.cxx
  vtkHierarchicalFractal.cxx
  vtkImageCompressor.cxx
  vtkIntegrateAttributes.cxx
  vtkIntegrateFlowThroughSurface.cxx
  vtkInteractorStyleTransferFunctionEditor.cxx
  vtkIntersectFragments.cxx
  vtkIsoVolume.cxx
  vtkKdTreeGenerator.cxx
  vtkKdTreeManager.cxx
  vtkMarkSelectedRows.cxx
  vtkMaterialInterfaceCommBuffer.cxx
  vtkMaterialInterfaceFilter.cxx
  vtkMaterialInterfaceIdList.cxx
  vtkMaterialInterfacePieceLoading.cxx
  vtkMaterialInterfacePieceTransaction.cxx
  vtkMaterialInterfacePieceTransactionMatrix.cxx
  vtkMaterialInterfaceProcessLoading.cxx
  vtkMaterialInterfaceProcessRing.cxx
  vtkMaterialInterfaceToProcMap.cxx
  vtkMergeArrays.cxx
  vtkMergeCompositeDataSet.cxx
  vtkMHTReader.cxx
  vtkMinMax.cxx
  vtkMPICompositeManager.cxx
  vtkMOCReader.cxx
  vtkMultiBlockTemporalStatistics.cxx
  vtkMultiProcessControllerHelper.cxx
  vtkOrderedCompositeDistributor.cxx
  vtkParallelSerialWriter.cxx
  vtkPConvertSelection.cxx
  vtkPEnSightGoldBinaryReader.cxx
  vtkPEnSightGoldReader.cxx
  vtkPEnSightReader.cxx
  vtkPExtractHistogram.cxx
  vtkPGenericEnSightReader.cxx
  vtkPhastaReader.cxx
  vtkPlotEdges.cxx
  vtkPointHandleRepresentationSphere.cxx
  vtkPolyLineToRectilinearGridFilter.cxx
  vtkPOPMatrices.cxx
  vtkPPhastaReader.cxx
  vtkPriorityHelper.cxx
  vtkPSciVizContingencyStats.cxx
  vtkPSciVizDescriptiveStats.cxx
  vtkPSciVizKMeans.cxx
  vtkPSciVizMultiCorrelativeStats.cxx
  vtkPSciVizPCAStats.cxx
  vtkPVAMRDualClip.cxx
  vtkPVAnimationCue.cxx
  vtkPVAnimationScene.cxx
  vtkPVArrayCalculator.cxx
  vtkPVArrowSource.cxx
  vtkPVAxesActor.cxx
  vtkPVAxesWidget.cxx
  vtkPVBooleanKeyFrame.cxx
  vtkPVBox.cxx
  vtkPVCenterAxesActor.cxx
  vtkPVClipClosedSurface.cxx
  vtkPVClipDataSet.cxx
  vtkPVCompositeDataPipeline.cxx
  vtkPVCompositeKeyFrame.cxx
  vtkPVConnectivityFilter.cxx
  vtkPVContourFilter.cxx
  vtkPVCueManipulator.cxx
  vtkPVDefaultPass.cxx
  vtkPVDReader.cxx
  vtkPVEnSightMasterServerReader.cxx
  vtkPVEnSightMasterServerReader2.cxx
  vtkPVEnSightMasterServerTranslator.cxx
  vtkPVExponentialKeyFrame.cxx
  vtkPVExtentTranslator.cxx
  vtkPVExtractSelection.cxx
  vtkPVExtractVOI.cxx
  vtkPVFrustumActor.cxx
  vtkPVGenericRenderWindowInteractor.cxx
  vtkPVGeometryFilter.cxx
  vtkPVGlyphFilter.cxx
  vtkPVHardwareSelector.cxx
  vtkPVInteractorStyle.cxx
  vtkPVJoystickFly.cxx
  vtkPVJoystickFlyIn.cxx
  vtkPVJoystickFlyOut.cxx
  vtkPVKeyFrame.cxx
  vtkPVKeyFrameAnimationCue.cxx
  vtkPVKeyFrameCueManipulator.cxx
  vtkPVLinearExtrusionFilter.cxx
  vtkPVLODActor.cxx
  vtkPVLODVolume.cxx
  vtkPVMergeTables.cxx
  vtkPVMultiBlockTemporalStatisticsExecutive.cxx
  vtkPVNullSource.cxx
  vtkPVPlane.cxx
  vtkPVPlotTime.cxx
  vtkPVPostFilter.cxx
  vtkPVPostFilterExecutive.cxx
  vtkPVRampKeyFrame.cxx
  vtkPVRecoverGeometryWireframe.cxx
  vtkPVRenderViewProxy.cxx
  vtkPVScalarBarActor.cxx
  vtkPVSelectionSource.cxx
  vtkPVSinusoidKeyFrame.cxx
  vtkPVTextSource.cxx
  vtkPVTrackballMoveActor.cxx
  vtkPVTrackballMultiRotate.cxx
  vtkPVTrackballPan.cxx
  vtkPVTrackballRoll.cxx
  vtkPVTrackballRotate.cxx
  vtkPVTrackballZoom.cxx
  vtkPVTransform.cxx
  vtkPVTrivialExtentTranslator.cxx
  vtkPVTrivialProducer.cxx
  vtkPVUpdateSuppressor.cxx
  vtkQuerySelectionSource.cxx
  vtkRealtimeAnimationPlayer.cxx
  vtkRectilinearGridConnectivity.cxx
  vtkReductionFilter.cxx
  vtkScatterPlotMapper.cxx
  vtkScatterPlotPainter.cxx
  vtkSciVizStatistics.cxx
  vtkSelectionConverter.cxx
  vtkSelectionSerializer.cxx
  vtkSequenceAnimationPlayer.cxx
  vtkSortedTableStreamer.cxx
  vtkSpyPlotBlock.cxx
  vtkSpyPlotBlockIterator.cxx
  vtkSpyPlotHistoryReader.cxx
  vtkSpyPlotIStream.cxx
  vtkSpyPlotReader.cxx
  vtkSpyPlotReaderMap.cxx
  vtkSpyPlotUniReader.cxx
  vtkSquirtCompressor.cxx
  vtkSurfaceVectors.cxx
  vtkTableFFT.cxx
  vtkTexturePainter.cxx
  vtkTilesHelper.cxx
  vtkTileDisplayHelper.cxx
  vtkTimestepsAnimationPlayer.cxx
  vtkTimeToTextConvertor.cxx
  vtkTrackballPan.cxx
  vtkTransferFunctionEditorRepresentation.cxx
  vtkTransferFunctionEditorRepresentation1D.cxx
  vtkTransferFunctionEditorRepresentationShapes1D.cxx
  vtkTransferFunctionEditorRepresentationShapes2D.cxx
  vtkTransferFunctionEditorRepresentationSimple1D.cxx
  vtkTransferFunctionEditorWidget.cxx
  vtkTransferFunctionEditorWidget1D.cxx
  vtkTransferFunctionEditorWidgetShapes1D.cxx
  vtkTransferFunctionEditorWidgetShapes2D.cxx
  vtkTransferFunctionEditorWidgetSimple1D.cxx
  vtkTransferFunctionViewer.cxx
  vtkUndoElement.cxx
  vtkUndoSet.cxx
  vtkUndoStack.cxx
  vtkUnstructuredPOPReader.cxx
  vtkUpdateSuppressorPipeline.cxx
  vtkVolumeRepresentationPreprocessor.cxx
  vtkVRMLSource.cxx
  vtkXMLCollectionReader.cxx
  vtkXMLPVDWriter.cxx
  vtkZlibImageCompressor.cxx
)

IF (VTK_USE_MPI)
  SET (Kit_SRCS
    ${Kit_SRCS}
    vtkAllToNRedistributeCompositePolyData.cxx
    vtkAllToNRedistributePolyData.cxx
    vtkBalancedRedistributePolyData.cxx
    vtkRedistributePolyData.cxx
    vtkWeightedRedistributePolyData.cxx
    )
  IF (PARAVIEW_USE_ICE_T)
    SET(Kit_SRCS
      ${Kit_SRCS}
      vtkIceTCompositePass.cxx
      vtkIceTContext.cxx)
    INCLUDE_DIRECTORIES(${ICET_INCLUDE_DIR})
    SET(KIT_LIBS ${KIT_LIBS} IceTCore IceTMPI IceTGL)
  ENDIF (PARAVIEW_USE_ICE_T)
  ADD_DEFINITIONS(-DPARAVIEW_USE_MPI)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (VTK_USE_MPI)

IF (VTK_USE_QVTK)
  INCLUDE_DIRECTORIES(
    ${QT_INCLUDE_DIR}
    ${ParaView_SOURCE_DIR}/VTK/GUISupport/Qt
    ${ParaView_BINARY_DIR}/VTK/GUISupport/Qt)
ENDIF (VTK_USE_QVTK)

IF(VTK_USE_SYSTEM_ZLIB)
  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
ENDIF(VTK_USE_SYSTEM_ZLIB)

SET_SOURCE_FILES_PROPERTIES(
  vtkAMRDualGridHelper.cxx
  vtkCacheSizeKeeper.cxx
  vtkMaterialInterfaceCommBuffer.cxx
  vtkMaterialInterfaceIdList.cxx
  vtkMaterialInterfacePieceLoading.cxx
  vtkMaterialInterfacePieceTransaction.cxx
  vtkMaterialInterfacePieceTransactionMatrix.cxx
  vtkMaterialInterfaceProcessLoading.cxx
  vtkMaterialInterfaceProcessRing.cxx
  vtkMaterialInterfaceToProcMap.cxx
  vtkPOPMatrices.cxx
  vtkPVPlotTime.cxx
  vtkSpyPlotBlock.cxx
  vtkSpyPlotBlockIterator.cxx
  vtkSpyPlotIStream.cxx
  vtkSpyPlotReaderMap.cxx
  vtkSpyPlotUniReader.cxx
  WRAP_EXCLUDE
  )

SET_SOURCE_FILES_PROPERTIES(
  vtkAbstractPOPReader.cxx
  vtkAnimationPlayer.cxx
  vtkCacheSizeKeeper.cxx
  vtkCameraManipulatorGUIHelper.cxx
  vtkCellIntegrator.cxx
  vtkImageCompressor.cxx
  vtkPEnSightReader.cxx
  vtkPVAnimationCue.cxx
  vtkPVCueManipulator.cxx
  vtkPVJoystickFly.cxx
  vtkPVKeyFrameAnimationCue.cxx
  vtkPVRenderViewProxy.cxx
  vtkSciVizStatistics.cxx
  vtkTileDisplayHelper.cxx
  vtkTransferFunctionEditorRepresentation.cxx
  vtkTransferFunctionEditorRepresentation1D.cxx
  vtkTransferFunctionEditorWidget.cxx
  vtkTransferFunctionEditorWidget1D.cxx
  vtkUndoElement.cxx

  ABSTRACT)

# The set of source files to be encoded.
SET(glsl_files
  vtkIceTCompositeZPassShader_fs
  )

INCLUDE_DIRECTORIES(
  ${PVCommon_SOURCE_DIR}
  ${PVCommon_BINARY_DIR}
  ${PVVTKExtensions_SOURCE_DIR}
  ${PVVTKExtensions_BINARY_DIR}
  ${VTK_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}
  ${SPCTH_INCLUDE_DIR}
  ${XDMF_INCLUDE_DIRS}
  )

# The following was copy/paste from vtkIncludeDirectories.cmake
# OpenGL include directories.
IF(APPLE)
  IF(VTK_USE_X OR VTK_OPENGL_HAS_OSMESA)
    INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIR})
  ENDIF(VTK_USE_X OR VTK_OPENGL_HAS_OSMESA)
ELSE(APPLE)
  INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIR})
ENDIF(APPLE)


SET(VTK_WRAP_HINTS ${CMAKE_CURRENT_SOURCE_DIR}/hints)

# Setup vtkInstantiator registration for this library's classes.
INCLUDE(${ParaView_SOURCE_DIR}/VTK/CMake/vtkMakeInstantiator.cmake)
VTK_MAKE_INSTANTIATOR3(vtk${KIT}Instantiator KitInstantiator_SRCS
  "${Kit_SRCS}"
  VTK_EXPORT
  "${CMAKE_CURRENT_BINARY_DIR}" "")

# load the ClientServer command
INCLUDE(${ParaView_SOURCE_DIR}/Utilities/VTKClientServer/vtkClientServer.cmake)
CS_INITIALIZE_WRAP()
INCLUDE_DIRECTORIES(${VTKCLIENTSERVER_INCLUDE_DIR})

# Wrap PVServer
# Create ClientServer wrappers for ParaView classes.
SET(_VTK_BINARY_DIR "${VTK_BINARY_DIR}")
SET(VTK_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")
INCLUDE(${VTK_CMAKE_DIR}/vtkExportKit.cmake)
VTK_EXPORT_KIT("${KIT}" "${KIT_UPCASE}" "${Kit_SRCS}")
SET(VTK_BINARY_DIR "${_VTK_BINARY_DIR}")

# Wrap PVServerCommon
SET(KIT_WRAP_DEPS PVCommon Parallel Hybrid Rendering VolumeRendering Widgets)
SET(VTK_KITS_DIR "${CMAKE_CURRENT_BINARY_DIR}/Utilities")
PV_WRAP_VTK_CS("${KIT}" "${KIT_UPCASE}" "${KIT_WRAP_DEPS}")

# Wrap for python
IF (PARAVIEW_ENABLE_PYTHON)
  SET(Kit_PYTHON_EXTRA_SRCS)
  SET(KIT_PYTHON_LIBS
    vtkPVCommonPythonD
    vtkHybridPythonD
    vtkWidgetsPythonD
    vtkParallelPythonD
    )
  # Tell vtkWrapPython.cmake to set VTK_PYTHON_LIBRARIES for us.
  SET(VTK_WRAP_PYTHON_FIND_LIBS 1)
  SET(VTK_WRAP_INCLUDE_DIRS ${VTK_INCLUDE_DIR})
  INCLUDE("${VTK_CMAKE_DIR}/vtkWrapPython.cmake")
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  SET(KIT_PYTHON_DEPS)

  SET(VTK_INSTALL_LIB_DIR "${VTKPythonWrapping_INSTALL_LIB_DIR}")
  SET(VTK_INSTALL_BIN_DIR "${VTKPythonWrapping_INSTALL_BIN_DIR}")
  SET(VTK_INSTALL_LIB_DIR_CM24 "${VTKPythonWrapping_INSTALL_LIB_DIR}")
  SET(VTK_INSTALL_BIN_DIR_CM24 "${VTKPythonWrapping_INSTALL_BIN_DIR}")

  # include python wrapping from VTK
  INCLUDE(KitCommonPythonWrapBlock)
ENDIF (PARAVIEW_ENABLE_PYTHON)


IF(VTK_USE_RENDERING)
  # Create custom commands to encode each assembly file into a C string
  # literal in a header file.

  SET(glsl_h_files)
  FOREACH(file ${glsl_files})
    SET(src ${CMAKE_CURRENT_SOURCE_DIR}/${file}.glsl)
    SET(res ${CMAKE_CURRENT_BINARY_DIR}/${file}.cxx)
    SET(resh ${CMAKE_CURRENT_BINARY_DIR}/${file}.h)
    SET(glsl_h_files ${glsl_h_files} "${resh}")
    ADD_CUSTOM_COMMAND(
      OUTPUT ${res} ${resh}
      DEPENDS ${src} vtkEncodeString
      COMMAND ${VTK_ENCODESTRING_EXE}
      ARGS ${res} ${src} ${file} --build-header VTK_EXPORT vtkSystemIncludes.h
      )
    SET(Kit_SRCS ${Kit_SRCS} ${res})
    SET_SOURCE_FILES_PROPERTIES( ${res} WRAP_EXCLUDE)
  ENDFOREACH(file)
ENDIF(VTK_USE_RENDERING)


# Create the library.
VTK_ADD_LIBRARY(vtk${KIT}
  ${Kit_SRCS}
  ${KitInstantiator_SRCS}
  )
TARGET_LINK_LIBRARIES(vtk${KIT}
  vtkPVCommon
  vtkPVCommonCS
  vtkXdmfCS
  vtkHybridCS
  vtkParallelCS
  vtkVolumeRenderingCS
  vtkWidgetsCS
  vtkViews
  vtkViewsCS
  vtkCharts
  vtkChartsCS
  KWCommon
  vtkNetCDF
  vtksys
  ${KIT_LIBS}
  ${PARAVIEW_HDF5_LIBRARIES}
  ${ZLIB_LIBRARY}
  ${SPCTH_LIBRARIES}
  ${OPENGL_gl_LIBRARY}
  ${PARAVIEW_ADDITIONAL_LIBRARIES}
  ${R_LIBRARIES}
)
IF (VTK_USE_QVTK)
  TARGET_LINK_LIBRARIES(vtk${KIT} QVTK)
ENDIF (VTK_USE_QVTK)

IF(NOT PV_INSTALL_NO_LIBRARIES)
  SET(temp_TO_INSTALL
    vtk${KIT})

  INSTALL(TARGETS ${temp_TO_INSTALL}
    EXPORT ${PV_INSTALL_EXPORT_NAME}
    RUNTIME DESTINATION ${PV_INSTALL_BIN_DIR} COMPONENT Runtime
    LIBRARY DESTINATION ${PV_INSTALL_LIB_DIR} COMPONENT Runtime
    ARCHIVE DESTINATION ${PV_INSTALL_LIB_DIR} COMPONENT Development)
ENDIF(NOT PV_INSTALL_NO_LIBRARIES)

# Development files install
IF (NOT PV_INSTALL_NO_DEVELOPMENT)
  GLOB_INSTALL_DEVELOPMENT(
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${PV_INSTALL_INCLUDE_DIR}
      "*.h;*.hxx;*.txx")
ENDIF (NOT PV_INSTALL_NO_DEVELOPMENT)

IF (BUILD_TESTING)
  SUBDIRS(Testing)
ENDIF (BUILD_TESTING)
