PROJECT(CorrespondenceEvaluator)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
    "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

SET( SOURCES ${SOURCES}
itkCorrespondenceEvaluator.h itkCorrespondenceEvaluator.txx)


ADD_EXECUTABLE(CorrespondenceEvaluator ${SOURCES} main.cpp)
TARGET_LINK_LIBRARIES(  CorrespondenceEvaluator
			ITKNumerics ITKCommon ITKBasicFilters ITKStatistics ITKIO
)

