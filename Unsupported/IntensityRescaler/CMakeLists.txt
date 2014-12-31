project(IntensityRescaler)

cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)




find_package(ITK REQUIRED)
if(ITK_FOUND)
   include(${ITK_USE_FILE})
else(ITK_FOUND)
   message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
endif(ITK_FOUND)

find_package(VTK REQUIRED)
if(VTK_FOUND)
  include(${VTK_USE_FILE})
else(VTK_FOUND)
  message(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
endif(VTK_FOUND)

include(CMake/qtSetup.cmake)
#find_package( QT REQUIRED )
#Find_Package( QT )
#if(QT_FOUND)
#  include(${QT_USE_FILE})
#else(QT_FOUND)
#  message(FATAL_ERROR, "Qt not found. Please set QT_DIR.")
#endif(QT_FOUND)


#Find_package( GLUT )
#Find_package( opengl )
include(CMake/glutSetup.cmake)

include(${CMAKE_ROOT}/Modules/FindOpenGL.cmake)
include(${CMAKE_ROOT}/Modules/FindGLU.cmake)


include_directories(
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

if(WIN32)
 #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /subsystem:windows")
endif(WIN32)

set(IntensityRescaler_MOC_SRCS
   IntensityRescalerGUIControls.h
)

if(QT_WRAP_CPP)
  qt_wrap_cpp(IntensityRescaler IntensityRescalerMoced_SRCS ${IntensityRescaler_MOC_SRCS})
endif(QT_WRAP_CPP)

set(IntensityRescaler_SRCS
  ${IntensityRescalerMoced_SRCS}
  IntensityRescaler.cxx
  BatchControls.cxx
  IntensityRescalerGUIControls.cxx
  IntensityRescalerCommandLine.cxx
)

set(IntensityRescaler_GUI_SRCS
   Batch.ui
   HistoGUI.ui
   IntensityRescalerGUI.ui
)




if(QT_WRAP_UI)
	qt_wrap_ui(IntensityRescaler IntensityRescaler_HDRS  IntensityRescaler_SRCS  ${IntensityRescaler_GUI_SRCS} )
endif(QT_WRAP_UI)



add_definitions(-DQT_DLL)
add_definitions(-DQT_THREAD_SUPPORT)

add_executable(IntensityRescaler ${IntensityRescaler_SRCS})

target_link_libraries(IntensityRescaler
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
