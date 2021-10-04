ttk_add_base_library(persistenceDiagramDictDecoding
  SOURCES
    PersistenceDiagramDictDecoding.cpp
  HEADERS
    PersistenceDiagramDictDecoding.h
  DEPENDS
    common
    persistenceDiagramClustering
  )
  target_compile_options(persistenceDiagramDictDecoding PRIVATE -O0 -ggdb
  )
