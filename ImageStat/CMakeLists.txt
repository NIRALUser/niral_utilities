CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(VERSION 2.8)

PROJECT(ImageStat)

FIND_PACKAGE(ITK REQUIRED)
IF (ITK_FOUND)
   INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
ENDIF (ITK_FOUND)

ADD_EXECUTABLE( ImageStat ImageStat.cc)
TARGET_LINK_LIBRARIES( ImageStat ${ITK_LIBRARIES} )

if( NOT ${ITK_USE_REVIEW_STATISTICS} )
  message(FATAL_ERROR "ITK has been compiled with ITK_USE_REVIEW_STATISTICS
set to OFF. ImageStat will fail to compile if this option is OFF. Please
recompile ITK with this option set to ON to compile ImageStat")
endif( NOT ${ITK_USE_REVIEW_STATISTICS} )

install( TARGETS ImageStat RUNTIME DESTINATION bin )
