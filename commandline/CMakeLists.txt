set(_SAVED_HDF5_DIR ${HDF5_DIR})
find_package(HDF5 1.8.5 NO_MODULE)
if(NOT HDF5_FOUND)
 set(HDF5_DIR ${_SAVED_HDF5_DIR} CACHE PATH "HDF5 install dir" FORCE)
 find_package(HDF5 1.8.5 REQUIRED)
endif()

find_package(Boost 1.54 REQUIRED COMPONENTS
  system filesystem program_options
)

add_executable(seamass-td
  seamass-td.cpp
)

target_include_directories(seamass-td PUBLIC
  ${HDF5_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
)

target_link_libraries(seamass-td
  seamass_topdown
  ${HDF5_LIBRARIES}
  ${Boost_LIBRARIES}
)

install(TARGETS seamass-td DESTINATION bin)
