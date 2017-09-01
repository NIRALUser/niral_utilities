project(TextureBioMarkerTool)

cmake_minimum_required(VERSION 3.5)

find_package(ITK COMPONENTS
  ITKAnisotropicSmoothing
  ITKBinaryMathematicalMorphology
  ITKCommon
  ITKConnectedComponents
  ITKIOGDCM
  ITKIOImageBase
  ITKIOTransformBase
  ITKImageFunction
  ITKImageGrid
  ITKImageIntensity
  ITKImageStatistics
  ITKMathematicalMorphology
  ITKRegionGrowing
  ITKReview
  ITKSmoothing
  ITKSpatialObjects
  ITKStatistics
  ITKThresholding
  ITKTransform
  ${ITK_IO_MODULES_USED}
)
add_executable(TextureBioMarkerTool TextureBioMarkerTool.cxx)
target_link_libraries(TextureBioMarkerTool ${ITK_LIBRARIES})

