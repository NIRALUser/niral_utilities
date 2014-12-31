project(DMDBioMarkerTool)
#find_path(LAPACK_LIBRARY_DIR liblapack.a)
cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

find_package(ITK REQUIRED)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
else(ITK_FOUND)
  message(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
endif(ITK_FOUND)

find_package(F2C REQUIRED)
if(F2C_FOUND)
  include(${F2C_USE_FILE})
else(F2C_FOUND)
  message(FATAL_ERROR "f2c.h not found. Please set f2c.h.")
endif(F2C_FOUND)

find_package(POLYFIT REQUIRED)
if(POLYFIT_FOUND)
  include(${POLYFIT_USE_FILE})
else(F2C_FOUND)
  message(FATAL_ERROR "polyfit not found. Please set polyfit.")
endif(POLYFIT_FOUND)

add_executable(DMDBioMarkerTool DMDBioMarkerTool.cxx)
target_link_libraries(DMDBioMarkerTool ${ITK_LIBRARIES} lapack)

