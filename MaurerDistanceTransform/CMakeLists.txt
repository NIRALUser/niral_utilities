project(MaurerDistanceTransform)

# Find SlicerExecutionModel
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})
 
# Find ITK
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
  ITKIONRRD
  ITKIONIFTI
  ITKIOTIFF
  ${ITK_IO_MODULES_USED}
  )
include(${ITK_USE_FILE})

set(MaurerDistanceTransform_SOURCE MaurerDistanceTransform.cxx)

SEMMacroBuildCLI(
  NAME MaurerDistanceTransform
  EXECUTABLE_ONLY
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  RUNTIME_OUTPUT_DIRECTORY ${NIRAL_UTILITIES_BINARY_DIR}/bin
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
  INSTALL_RUNTIME_DESTINATION ${INSTALL_RUNTIME_DESTINATION}
  INSTALL_LIBRARY_DESTINATION ${INSTALL_LIBRARY_DESTINATION}
  INSTALL_ARCHIVE_DESTINATION ${INSTALL_ARCHIVE_DESTINATION}
  )
