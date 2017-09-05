project(CorrespondenceEvaluator)

find_package(ITK)
if(ITK_FOUND)
    include(${ITK_USE_FILE})
else()
    message(FATAL_ERROR
    "ITK not found. Please set ITK_DIR.")
endif()

if( ${ITK_VERSION_MAJOR} VERSION_GREATER 3 )
  message( FATAL_ERROR "CorrespondenceEvaluator needs ITKv3 to compile" )
endif()

niral_add_executable(
  NAME CorrespondenceEvaluator
  SRC main.cpp
  ADDITIONAL_SRCS
    itkCorrespondenceEvaluator.h
    itkCorrespondenceEvaluator.txx
  TARGET_LIBRARIES
    ITKNumerics
    ITKCommon
    ITKBasicFilters
    ITKStatistics
    ITKIO
    lapack
  )
