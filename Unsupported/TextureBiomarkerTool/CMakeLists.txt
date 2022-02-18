project(TextureBioMarkerTool)

find_package(ITK 5.3 COMPONENTS
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

