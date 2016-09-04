set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
find_package(IntelOpenMP REQUIRED)

#SET(USE_MKL_64BIT_LIB On)
find_package(MKL REQUIRED)

add_library(matrix
  MatrixSparseMKL.hpp
  MatrixSparseMKL.cpp
)

target_include_directories(matrix PUBLIC
  ${MKL_INCLUDE_DIR}
)

target_compile_definitions(matrix
  PUBLIC ${MKL_COMPILER_DEFINITIONS}
)

target_link_libraries(matrix
  ${MKL_LIBRARIES}
  ${IntelOpenMP_LIBRARIES}
)

