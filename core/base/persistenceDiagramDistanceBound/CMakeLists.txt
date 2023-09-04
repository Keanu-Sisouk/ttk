ttk_add_base_library(persistenceDiagramDistanceBound
  SOURCES
    PersistenceDiagramDistanceBound.cpp
  HEADERS
    PersistenceDiagramDistanceBound.h
  DEPENDS
    triangulation
    geometry
    )

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(persistenceDiagramDistanceBound PUBLIC TTK_ENABLE_EIGEN)
  target_link_libraries(persistenceDiagramDistanceBound PUBLIC Eigen3::Eigen)
endif()