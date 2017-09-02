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

niral_add_executable(
  NAME TextureBioMarkerTool
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  )

