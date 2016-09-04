add_library(seamass_core
  Matrix.hpp
  DataMS.hpp
  DataMS.cpp
  Basis.hpp
  Basis.cpp
  OptimiserASRL.hpp
  OptimiserASRL.cpp
  BSpline.hpp
  BSpline.cpp
  BasisBSpline.hpp
  BasisBSpline.cpp
  BasisBSplineMZ.hpp
  BasisBSplineMZ.cpp
  BasisBSplineScale.hpp
  BasisBSplineScale.cpp
)

target_link_libraries(seamass_core
  matrix
)

