
set(SOURCES
  avtPlainTextFileFormat.C
  avtPlainTextFileFormat.h
  avtPlainTextOptions.C
  avtPlainTextOptions.h
  PlainTextPluginInfo.h
  PlainTextPluginInfo.C
  PlainTextEnginePluginInfo.C
  PlainTextMDServerPluginInfo.C
  PlainTextCommonPluginInfo.C
)

ADD_VISIT_READER(VisItPlainTextReader "1.0" 
  VISIT_READER_TYPE "STSD"
  VISIT_READER_NAME "avtPlainTextFileFormat"     
  VISIT_READER_USES_OPTIONS ON
  SERVER_SOURCES ${SOURCES}
  )
