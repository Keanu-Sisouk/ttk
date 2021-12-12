ttk_add_base_library(initDictPersistenceDiagram
  SOURCES
    InitDictPersistenceDiagram.cpp
  HEADERS
    InitDictPersistenceDiagram.h
  DEPENDS
    common
    persistenceDiagramAuction
)

target_compile_options(initDictPersistenceDiagram
  PRIVATE -ggdb
)
