PROJECT(TextureBioMarkerTool)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(VERSION 2.8)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)
list( FIND ITK_MODULES_ENABLED ITKReview FOUNDITKReview )
if( FOUNDITKReview EQUAL -1 )
  message(FATAL_ERROR "TextureBiomarkerTool: ITK needs to be compiled with Module_ITKReview set to ON")
endif()
ADD_EXECUTABLE(TextureBioMarkerTool TextureBioMarkerTool.cxx)
TARGET_LINK_LIBRARIES(TextureBioMarkerTool ${ITK_LIBRARIES}) 

