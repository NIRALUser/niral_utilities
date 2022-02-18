project(MultiAtlasSeg)

find_package(ITK 5.3 REQUIRED)
include(${ITK_USE_FILE})

list( FIND ITK_MODULES_ENABLED ITKReview position)
if( "${position}" EQUAL "-1" )
  message( FATAL_ERROR "MultiAtlasSeg needs ITKv4 compiled with ITKReview enabled" )
endif()

niral_add_executable(
  NAME MultiAtlasSeg
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  )
