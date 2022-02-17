project(TransformDeformationField)

if( ${ITK_VERSION_MAJOR} VERSION_GREATER 3 )
  message( FATAL_ERROR "TransformDeformationField needs ITKv3 to compile" )
endif()

find_package(ITK 3 REQUIRED)
include(${ITK_USE_FILE})

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

SEMMacroBuildCLI(
  NAME TransformDeformationField
  ${STATIC}
  TARGET_LIBRARIES
    ITKIO
    ITKAlgorithms
    ITKBasicFilters
  INCLUDE_DIRECTORIES
    dtiProcessFiles
  )
