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

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(persistenceDiagramDictEncoding PRIVATE TTK_ENABLE_EIGEN)
  target_link_libraries(persistenceDiagramDictEncoding PRIVATE Eigen3::Eigen)
endif()

