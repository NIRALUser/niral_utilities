project(ImageMath)

if(CMAKE_COMPILER_2005)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
endif()

find_package(ITK 5.3 REQUIRED COMPONENTS
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
  ITKImageNoise
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

file(GLOB LIBRARY_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.txx ${CMAKE_CURRENT_SOURCE_DIR}/*.hxx ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

set_target_properties(ImageMath PROPERTIES
  PUBLIC_HEADER "${LIBRARY_HEADERS}"
  )

install( TARGETS ImageMath 
  EXPORT niral_utilitiesTargets
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
  )
# Add all targets to the build-tree export set

set_property(GLOBAL APPEND PROPERTY niral_utilities_LIBRARIES ImageMath)
