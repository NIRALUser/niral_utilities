# NeuroLib_USE_QT()
# NeuroLib_USE_ITK()
# NeuroLib_USE_VTK()
# NeuroLib_USE_GLUT()

find_package(ITK REQUIRED)
if(ITK_FOUND)
   include(${ITK_USE_FILE})
else()
   message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
endif()

find_package(VTK REQUIRED)
if(VTK_FOUND)
  include(${VTK_USE_FILE})
else()
  message(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
endif()

include(${NIRAL_UTILITIES_SOURCE_DIR}/CMake/qtSetup.cmake)
include(${NIRAL_UTILITIES_SOURCE_DIR}/CMake/glutSetup.cmake)

include_directories(
    	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Datatypes
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Common
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Plot
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Window2D
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Window3D
	${QT_INCLUDE_DIR}
  	${GLUT_INCLUDE_DIR}
     )

link_libraries(
 ITKCommon
 ITKBasicFilters
 ${QT_QT_LIBRARY}
 ${GLUT_LIBRARIES}
)


set(Display_MOC_SRCS
  Window2D/QtWindow2D.h
  Plot/IntensityCurveBox.h
)

if(QT_WRAP_CPP)
  qt_wrap_cpp(Display Displaymoced_SRCS ${Display_MOC_SRCS})
endif()

        add_definitions(-DQT_DLL)
        add_definitions(-DQT_THREAD_SUPPORT)

set(Display_SRCS
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

add_library(Display ${Display_SRCS} )

