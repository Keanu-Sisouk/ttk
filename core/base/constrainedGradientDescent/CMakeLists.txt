ttk_add_base_library(constrainedGradientDescent
  SOURCES
    ConstrainedGradientDescent.cpp
  HEADERS
    ConstrainedGradientDescent.h
  DEPENDS
    common
)

target_compile_options(constrainedGradientDescent PRIVATE -O0 -ggdb -fsanitize=address)
target_link_options(constrainedGradientDescent PRIVATE -fsanitize=address)
