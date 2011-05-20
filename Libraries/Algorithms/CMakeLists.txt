

FIND_PACKAGE(ITK REQUIRED)
IF (ITK_FOUND)
   INCLUDE(${USE_ITK_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
ENDIF (ITK_FOUND)

FIND_PACKAGE(VTK REQUIRED)
IF (VTK_FOUND)
  INCLUDE (${VTK_USE_FILE})
ELSE(VTK_FOUND)
  MESSAGE(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
ENDIF (VTK_FOUND)

include(${NIRAL_UTILITIES_SOURCE_DIR}/CMake/qtSetup.cmake)

INCLUDE_DIRECTORIES(
    	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Datatypes
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Algorithms
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/IO/Common
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Common
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Plot
	${QT_INCLUDE_DIR}
	${ITK_SOURCE_DIR}/Utilities/zlib
     )

LINK_LIBRARIES(
 ${ITK_LIBRARIES}
 ${QT_QT_LIBRARY}
)


#SET(Display_MOC_SRCS 
#)

#IF(QT_WRAP_CPP)
#  QT_WRAP_CPP(Display Displaymoced_SRCS ${Display_MOC_SRCS})
#ENDIF(QT_WRAP_CPP)

ADD_DEFINITIONS(-DQT_DLL)
ADD_DEFINITIONS(-DQT_THREAD_SUPPORT)

SET(Algorithms_SRCS 
  ${Algorithms_SRCS}
  #${Algorithmsmoced_SRCS} 
  ImageHistogram.cxx
  ImageIntensityNormalizer.cxx
  ImageCombinaison.cxx
  ImageRotation.cxx
  Rotation3D.cxx
)

ADD_LIBRARY(Algorithms ${Algorithms_SRCS} )

