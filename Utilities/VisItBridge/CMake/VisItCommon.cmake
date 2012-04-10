
MACRO(VISIT_INSTALL_TARGETS target)     
  if(NOT PV_INSTALL_NO_LIBRARIES)    
    install(TARGETS ${target}
      EXPORT ${PV_INSTALL_EXPORT_NAME}
      RUNTIME DESTINATION ${PV_INSTALL_BIN_DIR} COMPONENT Runtime
      LIBRARY DESTINATION ${PV_INSTALL_LIB_DIR} COMPONENT Runtime
      ARCHIVE DESTINATION ${PV_INSTALL_LIB_DIR} COMPONENT Development)
  endif(NOT PV_INSTALL_NO_LIBRARIES)
  # Development
  if(NOT PV_INSTALL_NO_DEVELOPMENT)
   #VisIt has subdirectories so we need
   #to glob everything in the root directory
   #and all the headers in all the subdirectores
    GLOB_RECURSIVE_INSTALL_DEVELOPMENT(
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${PV_INSTALL_INCLUDE_DIR}
      "*.h;*.hxx;*.txx")
    if(dynamicHeaders)
      INSTALL(
        FILES ${dynamicHeaders}
        DESTINATION ${PV_INSTALL_INCLUDE_DIR}
        COMPONENT Development)
    endif(dynamicHeaders)
  endif(NOT PV_INSTALL_NO_DEVELOPMENT) 
  
ENDMACRO(VISIT_INSTALL_TARGETS)

FUNCTION(ADD_PARALLEL_LIBRARY target)
    VTK_ADD_LIBRARY(${target} ${ARGN})
    IF(VISIT_PARALLEL_CXXFLAGS)
        SET(PAR_COMPILE_FLAGS "")
        FOREACH (X ${VISIT_PARALLEL_CXXFLAGS})
            SET(PAR_COMPILE_FLAGS "${PAR_COMPILE_FLAGS} ${X}")
        ENDFOREACH(X)
        SET_TARGET_PROPERTIES(${target} PROPERTIES
            COMPILE_FLAGS ${PAR_COMPILE_FLAGS}
        )
        IF(VISIT_PARALLEL_LINKER_FLAGS)
            SET(PAR_LINK_FLAGS "")
            FOREACH(X ${VISIT_PARALLEL_LINKER_FLAGS})
                SET(PAR_LINK_FLAGS "${PAR_LINK_FLAGS} ${X}")
            ENDFOREACH(X)
            SET_TARGET_PROPERTIES(${target} PROPERTIES
                LINK_FLAGS ${PAR_LINK_FLAGS}
            )
        ENDIF(VISIT_PARALLEL_LINKER_FLAGS)

        IF(VISIT_PARALLEL_RPATH)
            SET(PAR_RPATHS "")
            FOREACH(X ${CMAKE_INSTALL_RPATH})
                SET(PAR_RPATHS "${PAR_RPATHS} ${X}")
            ENDFOREACH(X)
            FOREACH(X ${VISIT_PARALLEL_RPATH})
                SET(PAR_RPATHS "${PAR_RPATHS} ${X}")
            ENDFOREACH(X)
            SET_TARGET_PROPERTIES(${target} PROPERTIES
                INSTALL_RPATH ${PAR_RPATHS}
            )
        ENDIF(VISIT_PARALLEL_RPATH)

        IF(NOT VISIT_NOLINK_MPI_WITH_LIBRARIES)
            TARGET_LINK_LIBRARIES(${target} ${VISIT_PARALLEL_LIBS})
        ENDIF(NOT VISIT_NOLINK_MPI_WITH_LIBRARIES)
    ENDIF(VISIT_PARALLEL_CXXFLAGS)
ENDFUNCTION(ADD_PARALLEL_LIBRARY)

MACRO(VISIT_VTK_THIRD_PARTY_INCLUDE upper lower)
  if(VTK_USE_SYSTEM_${upper})
    if(${upper}_INCLUDE_DIR)
      include_directories(${${upper}_INCLUDE_DIR})
    endif(${upper}_INCLUDE_DIR)
  else(VTK_USE_SYSTEM_${upper})
    include_directories(
      ${VTK_BINARY_DIR}/Utilities/${lower}
      ${VTK_SOURCE_DIR}/Utilities/${lower}
    )
  endif(VTK_USE_SYSTEM_${upper})
ENDMACRO(VISIT_VTK_THIRD_PARTY_INCLUDE)

#called from readers that are being built into paraview
FUNCTION(ADD_VISIT_READER NAME VERSION)
  set(PLUGIN_NAME "vtk${NAME}")
  set(PLUGIN_VERSION "${VERSION}")
  set(ARG_VISIT_READER_NAME)
  set(ARG_VISIT_INCLUDE_NAME)
  set(ARG_VISIT_READER_TYPE)  
  set(ARG_VISIT_READER_USES_OPTIONS OFF)
  set(ARG_SERVER_SOURCES)  

  PV_PLUGIN_PARSE_ARGUMENTS(ARG 
    "VISIT_READER_NAME;VISIT_INCLUDE_NAME;VISIT_READER_TYPE;VISIT_READER_USES_OPTIONS;SERVER_SOURCES"
      "" ${ARGN} )    
  #check reader types
  string(REGEX MATCH "^[SM]T[SM]D$" VALID_READER_TYPE ${ARG_VISIT_READER_TYPE})

  if ( NOT VALID_READER_TYPE)
    MESSAGE(FATAL_ERROR "Invalid Reader Type. Valid Types are STSD, STMD, MTSD, MTMD")
  endif()

  #if the user hasn't defined an include name, we presume the reader name
  #is also the include name
  if(NOT ARG_VISIT_INCLUDE_NAME)
    set(ARG_VISIT_INCLUDE_NAME ${ARG_VISIT_READER_NAME})
  endif()
 
  if(ARG_VISIT_READER_USES_OPTIONS)
    #determine the name of the plugin info class by removing the 
    #avt from the start and the FileFormat from the back
    string(REGEX REPLACE "^avt|FileFormat$" "" TEMP_NAME ${ARG_VISIT_READER_NAME})
    set(ARG_VISIT_PLUGIN_INFO_HEADER ${TEMP_NAME}PluginInfo)
    set(ARG_VISIT_PLUGIN_INFO_CLASS ${TEMP_NAME}CommonPluginInfo)
  endif()

  set(XML_NAME ${NAME})
  set(LIBRARY_NAME "vtkVisItDatabases")
  #need to generate the VTK class wrapper
  string(SUBSTRING ${ARG_VISIT_READER_TYPE} 0 2 READER_WRAPPER_TYPE)
  
  configure_file(
    ${VISIT_CMAKE_DIR}/VisItExport.h.in
    ${VISIT_DATABASE_BINARY_DIR}/${PLUGIN_NAME}Export.h @ONLY)  
    
  configure_file(
      ${VISIT_CMAKE_DIR}/VisIt${READER_WRAPPER_TYPE}.h.in
      ${VISIT_DATABASE_BINARY_DIR}/${PLUGIN_NAME}.h @ONLY)  
  configure_file(
      ${VISIT_CMAKE_DIR}/VisIt${READER_WRAPPER_TYPE}.cxx.in
      ${VISIT_DATABASE_BINARY_DIR}/${PLUGIN_NAME}.cxx @ONLY)
  
  set(reader_sources
  ${VISIT_DATABASE_BINARY_DIR}/${PLUGIN_NAME}.cxx
  ${VISIT_DATABASE_BINARY_DIR}/${PLUGIN_NAME}.h
    )
    
  #fix up the arg_server_sources path for compilation
  set(ABS_SERVER_SOURCES "")
  foreach(SRC_FILENAME ${ARG_SERVER_SOURCES})
    list(APPEND ABS_SERVER_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILENAME}")
  endforeach()
    
  set(VISIT_SERVER_SOURCES ${VISIT_SERVER_SOURCES} ${reader_sources} 
    CACHE INTERNAL "vtk classes to wrap for client server")
  set(VISIT_DB_SOURCES ${VISIT_DB_SOURCES} ${ABS_SERVER_SOURCES}
    CACHE INTERNAL "visit sources for readers") 
  set(VISIT_DB_INC_DIRS ${VISIT_DB_INC_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}
    CACHE INTERNAL "include directories")  
  
ENDFUNCTION(ADD_VISIT_READER)

#called from readers that are being built into paraview
FUNCTION(ADD_VISIT_INTERFACE_READER NAME VERSION)

  set(INTERFACE_NAME "vtk${NAME}")
  set(INTERFACE_VERSION "${VERSION}")
  set(ARG_VISIT_READER_NAMES)
  set(ARG_VISIT_READER_TYPES)
  set(ARG_VISIT_READER_INCLUDES)
  set(ARG_VISIT_INTERFACE_CALL)
  set(ARG_VISIT_INTERFACE_FILE)
  set(ARG_VISIT_INTERFACE_EXEMPT_CLASSES)
  set(ARG_VISIT_READER_FILE_PATTERN)
  set(ARG_SERVER_SOURCES)
  PV_PLUGIN_PARSE_ARGUMENTS(ARG 
  "VISIT_READER_NAMES;VISIT_READER_TYPES;VISIT_READER_INCLUDES;VISIT_INTERFACE_CALL;VISIT_INTERFACE_FILE;VISIT_INTERFACE_EXEMPT_CLASSES;SERVER_SOURCES"
    "" ${ARGN} )   
    
  if (NOT ARG_VISIT_INTERFACE_CALL OR NOT ARG_VISIT_INTERFACE_FILE )
    MESSAGE(FATAL_ERROR "The macro file for the file interface needs to be defined.")
  endif(NOT ARG_VISIT_INTERFACE_CALL OR NOT ARG_VISIT_INTERFACE_FILE)

  #if the user hasn't defined an include name, we presume the reader name
  #is also the include name
  if(NOT ARG_VISIT_INCLUDE_NAME)
    set(ARG_VISIT_INCLUDE_NAME ${ARG_VISIT_READER_NAME})
  endif()
  
  set(LIBRARY_NAME "vtkVisItDatabases")
    
  list(LENGTH ARG_VISIT_READER_NAMES NUM_READERS)
  foreach( index RANGE ${NUM_READERS})
    if ( index LESS NUM_READERS )
      list(GET ARG_VISIT_READER_NAMES ${index} ARG_VISIT_READER_NAME)
      list(GET ARG_VISIT_READER_TYPES ${index} ARG_VISIT_READER_TYPE)
      list(GET ARG_VISIT_READER_INCLUDES ${index} ARG_VISIT_INCLUDE_NAME)
      
      #need to set up the vars needed by the configures
      string(REGEX REPLACE "^avt|FileFormat$" "" TEMP_NAME ${ARG_VISIT_READER_NAME})
      set(PLUGIN_NAME "vtkVisIt${TEMP_NAME}Reader")
      set(XML_NAME "VisItVisIt${TEMP_NAME}Reader")
    
      
      #need to generate the VTK class wrapper
      string(SUBSTRING ${ARG_VISIT_READER_TYPE} 0 2 READER_WRAPPER_TYPE)
      
      #determine if this file is exempt from the interface CanReadFile macro    
      list(FIND ARG_VISIT_INTERFACE_EXEMPT_CLASSES ${ARG_VISIT_READER_NAME} EXEMPT_READER)
      if ( EXEMPT_READER EQUAL -1 )
        set(VISIT_READER_USES_INTERFACE ON)
      else( EXEMPT_READER EQUAL -1 )
        set(VISIT_READER_USES_INTERFACE OFF)    
      endif( EXEMPT_READER EQUAL -1 ) 
  
      #we have to configure the macro file 
      configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_VISIT_INTERFACE_FILE}.h
        ${VISIT_DATABASE_BINARY_DIR}/${PLUGIN_NAME}${ARG_VISIT_INTERFACE_FILE}.h @ONLY)
        
      configure_file(
        ${VISIT_CMAKE_DIR}/VisItExport.h.in
        ${VISIT_DATABASE_BINARY_DIR}/${PLUGIN_NAME}Export.h @ONLY)  
          
      configure_file(
        ${VISIT_CMAKE_DIR}/VisIt${READER_WRAPPER_TYPE}.h.in
        ${VISIT_DATABASE_BINARY_DIR}/${PLUGIN_NAME}.h @ONLY)
        
      configure_file(
        ${VISIT_CMAKE_DIR}/VisIt${READER_WRAPPER_TYPE}.cxx.in
        ${VISIT_DATABASE_BINARY_DIR}/${PLUGIN_NAME}.cxx @ONLY)
        
      set(reader_sources
        ${VISIT_DATABASE_BINARY_DIR}/${PLUGIN_NAME}.cxx
        ${VISIT_DATABASE_BINARY_DIR}/${PLUGIN_NAME}.h)
    
    
      set(VISIT_SERVER_SOURCES ${VISIT_SERVER_SOURCES} ${reader_sources} 
        CACHE INTERNAL "vtk classes to wrap for client server")
        
    endif(index LESS NUM_READERS)
  endforeach(index)
  
  #fix up the arg_server_sources path for compilation
  set(ABS_SERVER_SOURCES "")
  foreach(SRC_FILENAME ${ARG_SERVER_SOURCES})
    list(APPEND ABS_SERVER_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILENAME}")
  endforeach()
        
  set(VISIT_DB_SOURCES ${VISIT_DB_SOURCES} ${ABS_SERVER_SOURCES}
    CACHE INTERNAL "visit sources for readers") 
    
  set(VISIT_DB_INC_DIRS ${VISIT_DB_INC_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}
    CACHE INTERNAL "include directories")
  

ENDFUNCTION(ADD_VISIT_INTERFACE_READER)

#Used for readers that are plugins for paraview
FUNCTION(ADD_VISIT_PLUGIN_READER NAME VERSION)
set(PLUGIN_NAME "vtk${NAME}")
set(PLUGIN_VERSION "${VERSION}")
set(ARG_VISIT_READER_NAME)
set(ARG_VISIT_INCLUDE_NAME)
set(ARG_VISIT_READER_TYPE)
set(ARG_VISIT_READER_FILE_PATTERN)
set(ARG_VISIT_READER_USES_OPTIONS OFF)
set(ARG_SERVER_SOURCES)

PV_PLUGIN_PARSE_ARGUMENTS(ARG 
  "VISIT_READER_NAME;VISIT_INCLUDE_NAME;VISIT_READER_TYPE;VISIT_READER_FILE_PATTERN;VISIT_READER_USES_OPTIONS;SERVER_SOURCES"
    "" ${ARGN} )    
#check reader types
string(REGEX MATCH "^[SM]T[SM]D$" VALID_READER_TYPE ${ARG_VISIT_READER_TYPE})

if ( NOT VALID_READER_TYPE)
  MESSAGE(FATAL_ERROR "Invalid Reader Type. Valid Types are STSD, STMD, MTSD, MTMD")
endif()  

#if the user hasn't defined an include name, we presume the reader name
#is also the include name
if(NOT ARG_VISIT_INCLUDE_NAME)
  set(ARG_VISIT_INCLUDE_NAME ${ARG_VISIT_READER_NAME})
endif()

MESSAGE(STATUS "Generating wrappings for ${PLUGIN_NAME}")
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${VISITBRIDGE_INCLUDE_DIRS}
  )  

if(ARG_VISIT_READER_USES_OPTIONS)
  #determine the name of the plugin info class by removing the 
  #avt from the start and the FileFormat from the back
  string(REGEX REPLACE "^avt|FileFormat$" "" TEMP_NAME ${ARG_VISIT_READER_NAME})
  set(ARG_VISIT_PLUGIN_INFO_HEADER ${TEMP_NAME}PluginInfo)
  set(ARG_VISIT_PLUGIN_INFO_CLASS ${TEMP_NAME}CommonPluginInfo)
endif()

set(XML_NAME ${NAME})
set(LIBRARY_NAME ${NAME})
#need to generate the VTK class wrapper
string(SUBSTRING ${ARG_VISIT_READER_TYPE} 0 2 READER_WRAPPER_TYPE)
configure_file(
    ${VISIT_CMAKE_DIR}/VisItExport.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}Export.h @ONLY)      
configure_file(
    ${VISIT_CMAKE_DIR}/VisIt${READER_WRAPPER_TYPE}.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}.h @ONLY)  
configure_file(
    ${VISIT_CMAKE_DIR}/VisIt${READER_WRAPPER_TYPE}.cxx.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}.cxx @ONLY)
  
#generate server manager xml file  
configure_file(
  ${VISIT_CMAKE_DIR}/VisItSM.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}SM.xml @ONLY)

#generate reader xml 
configure_file(
  ${VISIT_CMAKE_DIR}/VisItGUI.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}GUI.xml @ONLY)  
  
  
set(reader_sources
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}.cxx
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}.h
    )  
set(reader_server_xml
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}SM.xml
  )
set(reader_client_xml
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}GUI.xml
  )

#add the vtk classes to the argument list
set(PV_ARGS ${ARGN})
list(APPEND PV_ARGS "SERVER_MANAGER_SOURCES;${reader_sources}")

#now we need to add the XML info
list(APPEND PV_ARGS "SERVER_MANAGER_XML;${reader_server_xml}")
list(APPEND PV_ARGS "GUI_RESOURCE_FILES;${reader_client_xml}")

ADD_PARAVIEW_PLUGIN( ${NAME} ${VERSION} ${PV_ARGS} )
ENDFUNCTION(ADD_VISIT_PLUGIN_READER)

FUNCTION(ADD_VISIT_INTERFACE_PLUGIN_READER NAME VERSION)

set(INTERFACE_NAME "vtk${NAME}")
set(INTERFACE_VERSION "${VERSION}")
set(ARG_VISIT_READER_NAMES)
set(ARG_VISIT_READER_TYPES)
set(ARG_VISIT_READER_INCLUDES)
set(ARG_VISIT_INTERFACE_CALL)
set(ARG_VISIT_INTERFACE_FILE)
set(ARG_VISIT_INTERFACE_EXEMPT_CLASSES)
set(ARG_VISIT_READER_FILE_PATTERN)
set(ARG_SERVER_SOURCES)
    
PV_PLUGIN_PARSE_ARGUMENTS(ARG 
  "VISIT_READER_NAMES;VISIT_READER_TYPES;VISIT_READER_INCLUDES;VISIT_INTERFACE_CALL;VISIT_INTERFACE_FILE;VISIT_INTERFACE_EXEMPT_CLASSES;VISIT_READER_FILE_PATTERN;SERVER_SOURCES"
    "" ${ARGN} )   
    
if ( NOT ARG_VISIT_INTERFACE_CALL OR NOT ARG_VISIT_INTERFACE_FILE )
  MESSAGE(FATAL_ERROR "The macro file for the file interface needs to be defined.")
endif()

    
message(STATUS "Generating wrappings for ${INTERFACE_NAME}")    
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${VISITBRIDGE_INCLUDE_DIRS}
  )

#check reader types
set(INTERFACE_SOURCES)
set(INTERFACE_SMXML)
set(INTERFACE_GUIXML)
list(LENGTH ARG_VISIT_READER_NAMES NUM_READERS)
foreach( index RANGE ${NUM_READERS})
  if ( index LESS NUM_READERS )
    list(GET ARG_VISIT_READER_NAMES ${index} ARG_VISIT_READER_NAME)
    list(GET ARG_VISIT_READER_TYPES ${index} ARG_VISIT_READER_TYPE)
    list(GET ARG_VISIT_READER_INCLUDES ${index} ARG_VISIT_INCLUDE_NAME)
        
    #need to set up the vars needed by the configures
    string(REGEX REPLACE "^avt|FileFormat$" "" TEMP_NAME ${ARG_VISIT_READER_NAME})
    set(PLUGIN_NAME "vtk${TEMP_NAME}Reader")
    set(XML_NAME "VisIt${TEMP_NAME}Reader")
            
    #need to generate the VTK class wrapper
    string(SUBSTRING ${ARG_VISIT_READER_TYPE} 0 2 READER_WRAPPER_TYPE)
    
    #determine if this file is exempt from the interface CanReadFile macro    
    list(FIND ARG_VISIT_INTERFACE_EXEMPT_CLASSES ${ARG_VISIT_READER_NAME} EXEMPT_READER)
    if ( EXEMPT_READER EQUAL -1 )
      set(VISIT_READER_USES_INTERFACE ON)
    else()
      set(VISIT_READER_USES_INTERFACE OFF)    
    endif()    
    
    #we have to configure the macro file 
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_VISIT_INTERFACE_FILE}.h
        ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}${ARG_VISIT_INTERFACE_FILE}.h @ONLY)  
    
    #configure the declspec header
    configure_file(
        ${VISIT_CMAKE_DIR}/VisItExport.h.in
        ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}Export.h @ONLY)      
    
    #configure the header and implementation
    configure_file(
        ${VISIT_CMAKE_DIR}/VisIt${READER_WRAPPER_TYPE}.h.in
        ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}.h @ONLY)  
    configure_file(
        ${VISIT_CMAKE_DIR}/VisIt${READER_WRAPPER_TYPE}.cxx.in
        ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}.cxx @ONLY)
        
      
      
    #generate server manager xml file  
    configure_file(
      ${VISIT_CMAKE_DIR}/VisItSM.xml.in
      ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}SM.xml @ONLY)

    #generate reader xml 
    configure_file(
      ${VISIT_CMAKE_DIR}/VisItGUI.xml.in
      ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}GUI.xml @ONLY)  
  
    LIST(APPEND INTERFACE_SOURCES 
      ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}.cxx
      ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}.h
      )
    LIST(APPEND INTERFACE_SMXML 
      ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}SM.xml
      )
    LIST(APPEND INTERFACE_GUIXML 
      ${CMAKE_CURRENT_BINARY_DIR}/${PLUGIN_NAME}GUI.xml
      )
    
  endif()
endforeach( index )    

#add the vtk classes to the argument list
set(PV_ARGS ${ARGN})
list(APPEND PV_ARGS "SERVER_MANAGER_SOURCES;${INTERFACE_SOURCES}")

#now we need to add the XML info
list(APPEND PV_ARGS "SERVER_MANAGER_XML;${INTERFACE_SMXML}")
list(APPEND PV_ARGS "GUI_RESOURCE_FILES;${INTERFACE_GUIXML}")

ADD_PARAVIEW_PLUGIN( ${NAME} ${VERSION} ${PV_ARGS} )

ENDFUNCTION(ADD_VISIT_INTERFACE_PLUGIN_READER)