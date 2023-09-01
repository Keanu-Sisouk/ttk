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
  target_compile_definitions(constrainedGradientDescent PRIVATE TTK_ENABLE_EIGEN)
  target_link_libraries(constrainedGradientDescent PRIVATE Eigen3::Eigen)
endif()