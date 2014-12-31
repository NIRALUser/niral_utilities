project(CorrespondenceEvaluator)

find_package(ITK)
if(ITK_FOUND)
    include(${ITK_USE_FILE})
else(ITK_FOUND)
    message(FATAL_ERROR
    "ITK not found. Please set ITK_DIR.")
endif(ITK_FOUND)

if( ${ITK_VERSION_MAJOR} VERSION_GREATER 3 )
  message( FATAL_ERROR "CorrespondenceEvaluator needs ITKv3 to compile" )
endif()


set( SOURCES ${SOURCES}
itkCorrespondenceEvaluator.h itkCorrespondenceEvaluator.txx)


add_executable(CorrespondenceEvaluator ${SOURCES} main.cpp)
target_link_libraries(  CorrespondenceEvaluator
			ITKNumerics ITKCommon ITKBasicFilters ITKStatistics ITKIO
)

