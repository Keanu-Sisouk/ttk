ttk_add_base_library(persistenceDiagramDictEncoding
  SOURCES
    PersistenceDiagramDictEncoding.cpp
    InitDictPersistenceDiagram.cpp
  HEADERS
    PersistenceDiagramDictEncoding.h
    InitDictPersistenceDiagram.h
  DEPENDS
    common
    persistenceDiagramAuction
    constrainedGradientDescent
    persistenceDiagramClustering
  )
  target_compile_options(persistenceDiagramDictEncoding PRIVATE -O0 -ggdb
  )
