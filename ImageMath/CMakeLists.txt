project(ImageMath)

cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

if(CMAKE_COMPILER_2005)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
endif()

find_package(ITK COMPONENTS
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
  ${ITK_IO_MODULES_USED}
)
include(${ITK_USE_FILE})

if( ${ITK_VERSION} VERSION_LESS 4.4.0 )
  message( FATAL_ERROR "ImageMath needs ITKv4.4.0 or greater to compile" )
endif()

add_executable(ImageMath ImageMath.cxx )
target_link_libraries( ImageMath ${ITK_LIBRARIES} )

install( TARGETS ImageMath RUNTIME DESTINATION bin )
