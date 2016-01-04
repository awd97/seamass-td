set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(MKL REQUIRED)

# gcc doesn't like SHARED build atm
#add_library(seamass_topdown SHARED
add_library(seamass_topdown STATIC
  topdown.hpp
  topdown.cpp
  BasisChargeDistribution.hpp
  BasisChargeDistribution.cpp
)

include(GenerateExportHeader)
generate_export_header(seamass_topdown)

target_include_directories(seamass_topdown PUBLIC
  ${CMAKE_BINARY_DIR}/topdown
)

target_link_libraries(seamass_topdown
  seamass_core
)

#install(TARGETS seamass_topdown
#        RUNTIME DESTINATION bin
# 		LIBRARY DESTINATION lib)
#install(FILES topdown.hpp DESTINATION seamass/include)
