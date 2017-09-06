project(MultiAtlasSeg)

cmake_minimum_required(VERSION 3.5)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

list( FIND ITK_MODULES_ENABLED ITKReview position)
if( "${position}" EQUAL "-1" )
  message( FATAL_ERROR "MultiAtlasSeg needs ITKv4 compiled with ITKReview enabled" )
endif()

niral_add_executable(
  NAME MultiAtlasSeg
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  )
