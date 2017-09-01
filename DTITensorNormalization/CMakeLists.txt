cmake_minimum_required(VERSION 3.5)
PROJECT(DTITensorNormalization)




FIND_PACKAGE(ITK)

IF (ITK_FOUND)
   include(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
ENDIF (ITK_FOUND)

ADD_EXECUTABLE( DTITensorNormalization DTITensorNormalization.cxx)


TARGET_LINK_LIBRARIES(DTITensorNormalization ${ITK_LIBRARIES} )

install( TARGETS DTITensorNormalization RUNTIME DESTINATION bin )

