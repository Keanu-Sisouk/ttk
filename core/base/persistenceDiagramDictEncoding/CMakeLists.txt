ttk_add_base_library(persistenceDiagramDictEncoding
  SOURCES
    PersistenceDiagramDictEncoding.cpp
    InitDictBorder.cpp
    InitDictRandomly.cpp
  HEADERS
    PersistenceDiagramDictEncoding.h
    InitDictBorder.h
    InitDictRandomly.h
  DEPENDS
    common
    persistenceDiagramAuction
    constrainedGradientDescent
    persistenceDiagramClustering
  )
  target_compile_options(persistenceDiagramDictEncoding PRIVATE -O0 -ggdb
  )
