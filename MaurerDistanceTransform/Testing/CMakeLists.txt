cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

project(MaurerDistanceTransform)

# Find SlicerExecutionModel
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})
 
# Find ITK
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
ITKTestKernel
ITKIONRRD
ITKIONIFTI
ITKIOTIFF
${ITK_IO_MODULES_USED}
)
include(${ITK_USE_FILE})

if( ${ITK_VERSION} VERSION_LESS 4.4.0 )
message( FATAL_ERROR "MaurerDistanceTransform needs ITKv4.4.0 or greater to compile" )
endif()


set(MaurerDistanceTransform_SOURCE MaurerDistanceTransform.cxx)

SEMMacroBuildCLI(
  NAME MaurerDistanceTransform
  EXECUTABLE_ONLY
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  )

