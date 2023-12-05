ttk_add_base_library(persistenceDiagramSlicedWasserstein
  SOURCES
    PersistenceDiagramSlicedWasserstein.cpp
  HEADERS
    PersistenceDiagramSlicedWasserstein.h
  DEPENDS
    common
    persistenceDiagramAuction
  )

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(persistenceDiagramSlicedWasserstein PRIVATE TTK_ENABLE_EIGEN)
  target_link_libraries(persistenceDiagramSlicedWasserstein PRIVATE Eigen3::Eigen)
endif()

