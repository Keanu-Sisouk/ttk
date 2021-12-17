ttk_add_base_library(constrainedGradientDescent
  SOURCES
    ConstrainedGradientDescent.cpp
  HEADERS
    ConstrainedGradientDescent.h
  DEPENDS
    common
)

target_compile_options(constrainedGradientDescent
  PRIVATE -ggdb
)
