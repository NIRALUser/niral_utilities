project(IntensityRescaler)

#snippet taken from: http://www.guyrutenberg.com/2014/01/05/enabling-c11-c0x-in-cmake/
include(CheckCXXCompilerFlag)
if( EXTRA_FLAG )
  set( CMAKE_CXX_FLAGS_LIST ${CMAKE_CXX_FLAGS} )
  separate_arguments(CMAKE_CXX_FLAGS_LIST)
  get_property( helpstring CACHE CMAKE_CXX_FLAGS PROPERTY HELPSTRING )
  list(FIND CMAKE_CXX_FLAGS_LIST ${EXTRA_FLAG} found_flag )
  if( found_flag LESS 0 )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAG}" CACHE STRING "${helpstring}" FORCE )
  endif()
endif()

find_package(ITK 5.3 REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK 9 REQUIRED)
include(${VTK_USE_FILE})


include_directories(
  Libraries/Datatypes
  Libraries/IO/Common
  Libraries/Algorithms
  Libraries/Display/Common
  Libraries/Display/Plot
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

option(COMPILE_IntensityRescaler_GUI "Qt IntensityRescaler GUI support." OFF)

if( COMPILE_IntensityRescaler_GUI )
  find_package( Qt REQUIRED )
  include(${QT_USE_FILE})

  include_directories(
    ${QT_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}
    Libraries/IO/Qt
    Libraries/Display/Window2D
    )
  if(WIN32)
    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /subsystem:windows")
  endif()
  #Find_package( GLUT )
  #Find_package( opengl )
  include(CMake/glutSetup.cmake)

  include(${CMAKE_ROOT}/Modules/FindOpenGL.cmake)
  include(${CMAKE_ROOT}/Modules/FindGLU.cmake)

  set(IntensityRescaler_MOC_SRCS
    IntensityRescalerGUIControls.h
    )
  if(QT_WRAP_CPP)
    qt_wrap_cpp(IntensityRescaler IntensityRescalerMoced_SRCS ${IntensityRescaler_MOC_SRCS})
  endif()
  set( IntensityRescaler_GUI_CXX_SRCS 
    IntensityRescalerGUIControls.cxx
    BatchControls.cxx
    )
  set(IntensityRescaler_GUI_SRCS
    Batch.ui
    HistoGUI.ui
    IntensityRescalerGUI.ui
    )
  if(QT_WRAP_UI)
	  qt_wrap_ui(IntensityRescaler IntensityRescaler_HDRS  IntensityRescaler_SRCS  ${IntensityRescaler_GUI_SRCS} )
  endif()
  add_definitions(-DQT_DLL)
  add_definitions(-DQT_THREAD_SUPPORT)
  add_definitions( -DQT_GUI )

  set(OPEN_GL_LIBRARIES
    ${QT_GL_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    ${OPENGL_LIBRARY}
    ${GLUT_LIBRARIES}
    )
  set( QT_LIBRARIES_LIST
    ${QT_QT_LIBRARY}
    ${QT_LIBRARIES}
    ${QT_QTMAIN_LIBRARY}
    )
else()
  add_definitions( -DQT_NO_GUI )
endif()
set(IntensityRescaler_SRCS
  ${IntensityRescalerMoced_SRCS}
  ${IntensityRescaler_GUI_CXX_SRCS}
  IntensityRescaler.cxx
  IntensityRescalerCommandLine.cxx
  )

add_subdirectory(Libraries)

add_executable(IntensityRescaler ${IntensityRescaler_SRCS})

target_link_libraries(IntensityRescaler
  Algorithms
  Display
  ${ITK_LIBRARIES}
  ${QT_LIBRARIES_LIST}
  ${OPEN_GL_LIBRARIES}
  ${VTK_LIBRARIES}
  )

install( TARGETS IntensityRescaler RUNTIME DESTINATION bin )

