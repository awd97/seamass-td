set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

SET(USE_MKL_64BIT_LIB On)
ADD_DEFINITIONS(-DMKL_ILP64)
find_package(MKL REQUIRED)

add_library(seamass_core
  SparseMatrixMKL.hpp
  SparseMatrixMKL.cpp
  core.hpp
  core.cpp
  BSpline.hpp
  BSpline.cpp
  Basis.hpp
  Basis.cpp
  OptimiserASRL.hpp
  OptimiserASRL.cpp
)

target_include_directories(seamass_core PUBLIC
  ${MKL_INCLUDE_DIR}
)

target_link_libraries(seamass_core
  ${MKL_LIBRARIES}
)
