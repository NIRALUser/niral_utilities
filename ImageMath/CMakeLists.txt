project(ImageMath)

cmake_minimum_required(VERSION 3.5)

if(CMAKE_COMPILER_2005)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
endif()

find_package(ITK 4.4.0 REQUIRED COMPONENTS
  ITKAnisotropicSmoothing
  ITKBinaryMathematicalMorphology
  ITKCommon
  ITKConnectedComponents
  ITKCurvatureFlow
  ITKDistanceMap
  ITKIOImageBase
  ITKImageFilterBase
  ITKImageFunction
  ITKImageGrid
  ITKImageIntensity
  ITKMathematicalMorphology
  ITKSmoothing
  ITKStatistics
  ITKThresholding
  ITKTransform
  ITKTestKernel
  ${ITK_IO_MODULES_USED}
  )
include(${ITK_USE_FILE})

add_library( ImageMathLib ImageMath.cxx )
set_target_properties( ImageMathLib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint" ) 
target_link_libraries( ImageMathLib ${ITK_LIBRARIES} )

niral_add_executable(
  NAME ImageMath
  SRC ImageMathApp.cxx
  ADDITIONAL_SRCS
    ImageMathApp.cxx
  TARGET_LIBRARIES
    ImageMathLib
  )

if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()
