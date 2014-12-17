PROJECT(IntensityRescaler)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(VERSION 2.8)




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

include(CMake/qtSetup.cmake)
#find_package( QT REQUIRED )
#Find_Package( QT )
#IF (QT_FOUND)
#  INCLUDE (${QT_USE_FILE})
#ELSE(QT_FOUND)
#  MESSAGE(FATAL_ERROR, "Qt not found. Please set QT_DIR.")
#ENDIF (QT_FOUND)


#Find_package( GLUT )
#Find_package( opengl )
include(CMake/glutSetup.cmake)

INCLUDE (${CMAKE_ROOT}/Modules/FindOpenGL.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/FindGLU.cmake)


INCLUDE_DIRECTORIES(
  ${QT_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  Libraries/Display/Plot
  Libraries/Datatypes
  Libraries/IO/Qt
  Libraries/IO/Common
  Libraries/Algorithms
  Libraries/Display/Common
  Libraries/Display/Window2D
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

IF(WIN32)
 #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /subsystem:windows")
ENDIF(WIN32)

SET(IntensityRescaler_MOC_SRCS 
   IntensityRescalerGUIControls.h
)

IF(QT_WRAP_CPP)
  QT_WRAP_CPP(IntensityRescaler IntensityRescalerMoced_SRCS ${IntensityRescaler_MOC_SRCS})
ENDIF(QT_WRAP_CPP)

SET(IntensityRescaler_SRCS 
  ${IntensityRescalerMoced_SRCS} 
  IntensityRescaler.cxx
  BatchControls.cxx
  IntensityRescalerGUIControls.cxx
  IntensityRescalerCommandLine.cxx
)

SET(IntensityRescaler_GUI_SRCS 
   Batch.ui
   HistoGUI.ui
   IntensityRescalerGUI.ui
)




IF(QT_WRAP_UI)
	QT_WRAP_UI(IntensityRescaler IntensityRescaler_HDRS  IntensityRescaler_SRCS  ${IntensityRescaler_GUI_SRCS} )
ENDIF(QT_WRAP_UI)



ADD_DEFINITIONS(-DQT_DLL)
ADD_DEFINITIONS(-DQT_THREAD_SUPPORT)

ADD_EXECUTABLE(IntensityRescaler ${IntensityRescaler_SRCS})

TARGET_LINK_LIBRARIES(IntensityRescaler
  Display
  Algorithms
  ${ITK_LIBRARIES}
  ${QT_QT_LIBRARY}
  ${QT_QTMAIN_LIBRARY}
  ${QT_GL_LIBRARY}
  ${OPENGL_glu_LIBRARY}
  ${OPENGL_LIBRARY}
  ${GLUT_LIBRARIES}
  vtkRendering
  vtkGraphics
  vtkImaging
  vtkIO
  vtkFiltering
  vtkCommon
  vtkHybrid
  #vtkPatented
 )
