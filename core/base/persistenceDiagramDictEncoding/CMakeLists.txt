ttk_add_base_library(persistenceDiagramDictEncoding
  SOURCES
    PersistenceDiagramDictEncoding.cpp
  HEADERS
    PersistenceDiagramDictEncoding.h
  DEPENDS
    common
    persistenceDiagramAuction
    constrainedGradientDescent
    persistenceDiagramClustering
    initDictPersistenceDiagram
  )

target_compile_options(persistenceDiagramDictEncoding PRIVATE -O0 -ggdb -fsanitize=address)
target_link_options(persistenceDiagramDictEncoding PRIVATE -fsanitize=address)
