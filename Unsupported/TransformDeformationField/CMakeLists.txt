
cmake_minimum_required(VERSION 3.5)

project(TransformDeformationField)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if( ${ITK_VERSION_MAJOR} VERSION_GREATER 3 )
  message( FATAL_ERROR "TransformDeformationField needs ITKv3 to compile" )
endif()

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

