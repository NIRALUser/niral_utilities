PROJECT(convertITKformats)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(VERSION 2.8)

FIND_PACKAGE(ITK REQUIRED)
IF (ITK_FOUND)
   INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
ENDIF (ITK_FOUND)

if( ${ITK_VERSION_MAJOR} VERSION_GREATER 3)
  if( NOT ITKV3_COMPATIBILITY )
    message(FATAL_ERROR "convertITKformats: ITK needs to be compiled with ITKV3_COMPATIBILITY set to ON")
  endif()
endif()

ADD_EXECUTABLE( convertITKformats convertITKformats.cxx)
ADD_EXECUTABLE( convertITKformats4D convertITKformats4D.cxx)
ADD_EXECUTABLE( convertTextToImage convertTextToImage.cxx)

TARGET_LINK_LIBRARIES( convertITKformats ${ITK_LIBRARIES} )
TARGET_LINK_LIBRARIES( convertITKformats4D ${ITK_LIBRARIES} )
TARGET_LINK_LIBRARIES( convertTextToImage ${ITK_LIBRARIES} )

install( TARGETS convertITKformats convertITKformats4D convertTextToImage RUNTIME DESTINATION bin )
