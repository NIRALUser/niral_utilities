# NeuroLib_USE_QT()
# NeuroLib_USE_ITK()
# NeuroLib_USE_VTK()
# NeuroLib_USE_GLUT()

FIND_PACKAGE(ITK REQUIRED)
IF (ITK_FOUND)
   INCLUDE(${ITK_USE_FILE})
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
include(${NIRAL_UTILITIES_SOURCE_DIR}/CMake/glutSetup.cmake)

INCLUDE_DIRECTORIES(
    	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Datatypes
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Common
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Plot
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Window2D
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Window3D
	${QT_INCLUDE_DIR}
  	${GLUT_INCLUDE_DIR}
     )

LINK_LIBRARIES(
 ITKCommon
 ITKBasicFilters
 ${QT_QT_LIBRARY}
 ${GLUT_LIBRARIES}
)


SET(Display_MOC_SRCS 
  Window2D/QtWindow2D.h
  Plot/IntensityCurveBox.h
)

IF(QT_WRAP_CPP)
  QT_WRAP_CPP(Display Displaymoced_SRCS ${Display_MOC_SRCS})
ENDIF(QT_WRAP_CPP)

        ADD_DEFINITIONS(-DQT_DLL)
        ADD_DEFINITIONS(-DQT_THREAD_SUPPORT)

SET(Display_SRCS 
  ${Display_SRCS}
  ${Displaymoced_SRCS} 
  Window2D/QtWindow2D.cxx
  Plot/FunctionPlot2D.cxx
  Plot/FunctionPlot2DBox.cxx
  Plot/BarGraph2D.cxx
  Plot/BarGraph2DBox.cxx
  Plot/IntensityCurveVTK.cxx
  Plot/IntensityCurveBox.cxx
#  Window3D/VolumeRenderingActor.cxx
#  Window3D/VolumeRenderingActor.h
#  Window3D/IsoContourActor.cxx
#  Window3D/IsoContourActor.h
)

ADD_LIBRARY(Display ${Display_SRCS} )

