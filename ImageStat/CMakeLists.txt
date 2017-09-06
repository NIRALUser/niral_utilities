cmake_minimum_required(VERSION 3.5)

project(ImageStat)

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

niral_add_executable(
  NAME ImageStat
  SRC ImageStat.cc
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  )
