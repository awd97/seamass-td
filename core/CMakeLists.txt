set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(MKL REQUIRED)

add_library(seamass_core STATIC
  core.hpp
  core.cpp
  Basis.hpp
  Basis.cpp
  BSpline.hpp
  BSpline.cpp
  OptimiserASRL.hpp
  OptimiserASRL.cpp
)

target_include_directories(seamass_core PUBLIC
  ${MKL_INCLUDE_DIR}
)

target_link_libraries(seamass_core
  ${MKL_LIBRARIES}
)