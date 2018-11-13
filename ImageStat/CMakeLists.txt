project(ImageStat)

cmake_minimum_required(VERSION 3.5)

find_package(ITK REQUIRED COMPONENTS
  ITKCommon
  ITKDistanceMap
  ITKIOImageBase
  ITKImageFilterBase
  ITKImageIntensity
  ITKStatistics
  ITKThresholding
  ${ITK_IO_MODULES_USED}
  )
include(${ITK_USE_FILE})

if( NOT ${ITK_USE_REVIEW_STATISTICS} )
  message(FATAL_ERROR "ITK has been compiled with ITK_USE_REVIEW_STATISTICS
set to OFF. ImageStat will fail to compile if this option is OFF. Please
recompile ITK with this option set to ON to compile ImageStat")
endif()

add_library( ImageStatLib ImageStat.cc )
target_link_libraries( ImageStatLib ${ITK_LIBRARIES} )


niral_add_executable(
  NAME ImageStat
  SRC ImageStat.cc
  TARGET_LIBRARIES
    ImageStatLib
  )

file(GLOB LIBRARY_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.txx ${CMAKE_CURRENT_SOURCE_DIR}/*.hxx ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

set_target_properties(ImageStat PROPERTIES
  PUBLIC_HEADER "${LIBRARY_HEADERS}"
  )


install( TARGETS ImageStat 
  EXPORT niral_utilitiesTargets
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
  )

set_property(GLOBAL APPEND PROPERTY niral_utilities_LIBRARIES ImageStat)