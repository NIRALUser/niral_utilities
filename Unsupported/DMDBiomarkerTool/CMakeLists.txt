project(DMDBioMarkerTool)
#find_path(LAPACK_LIBRARY_DIR liblapack.a)
cmake_minimum_required(VERSION 3.5)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if( ${ITK_VERSION_MAJOR} VERSION_GREATER 3 )
  message( FATAL_ERROR "TransformDeformationField needs ITKv3 to compile" )
endif()

find_package(F2C REQUIRED)
if(F2C_FOUND)
  include(${F2C_USE_FILE})
else()
  message(FATAL_ERROR "f2c.h not found. Please set f2c.h.")
endif()

find_package(POLYFIT REQUIRED)
if(POLYFIT_FOUND)
  include(${POLYFIT_USE_FILE})
else()
  message(FATAL_ERROR "polyfit not found. Please set polyfit.")
endif()

add_executable(DMDBioMarkerTool DMDBioMarkerTool.cxx)
target_link_libraries(DMDBioMarkerTool ${ITK_LIBRARIES} lapack)

