
if( COMPILE_IntensityRescaler_GUI )
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
 ${QT_LIBRARIES}
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

else()

include_directories(
    	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Datatypes
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Common
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Plot
     )
link_libraries(
   ${VTK_LIBRARIES}
vtkFiltering
)

set(Display_SRCS
  Plot/IntensityCurveVTK.cxx
)

add_library(Display ${Display_SRCS} )

endif()
