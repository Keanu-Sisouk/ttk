ttk_add_base_library(persistenceDiagramSlicedWasserstein
  SOURCES
    PersistenceDiagramSlicedWasserstein.cpp
  HEADERS
    PersistenceDiagramSlicedWasserstein.h
  DEPENDS
    common
    persistenceDiagramAuction
  )

target_compile_options(persistenceDiagramSlicedWasserstein PRIVATE -O0 -ggdb -fsanitize=address)
target_link_options(persistenceDiagramSlicedWasserstein PRIVATE -fsanitize=address)