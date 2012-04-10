
set(SOURCES
  avtTecplotFileFormat.C  
  TecplotFile.C
  TecplotFile.h
  )
  
set(BINARY_SOURCES
  avtTecplotBinaryFileFormat.C
  TecplotFile.C
  )
ADD_VISIT_READER(VisItTecplotReader "1.0"
  VISIT_READER_TYPE "STMD"
  VISIT_READER_NAME "avtTecplotFileFormat"    
  SERVER_SOURCES ${SOURCES}
  )
ADD_VISIT_READER(VisItTecplotBinaryReader "1.0"
  VISIT_READER_TYPE "STMD"
  VISIT_READER_NAME "avtTecplotBinaryFileFormat"    
  SERVER_SOURCES ${BINARY_SOURCES}
  )