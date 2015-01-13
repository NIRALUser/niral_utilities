

include_directories(
    	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Datatypes
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Algorithms
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/IO/Common
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Common
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Plot
	${QT_INCLUDE_DIR}
	${ITK_SOURCE_DIR}/Utilities/zlib
     )

link_libraries(
 ${ITK_LIBRARIES}
 ${QT_QT_LIBRARY}
 ${QT_LIBRARIES}
)


#set(Display_MOC_SRCS
#)

#if(QT_WRAP_CPP)
#  qt_wrap_cpp(Display Displaymoced_SRCS ${Display_MOC_SRCS})
#endif()

add_definitions(-DQT_DLL)
add_definitions(-DQT_THREAD_SUPPORT)

set(Algorithms_SRCS
  ${Algorithms_SRCS}
  #${Algorithmsmoced_SRCS}
  ImageHistogram.cxx
  ImageIntensityNormalizer.cxx
  ImageCombinaison.cxx
  ImageRotation.cxx
  Rotation3D.cxx
)

add_library(Algorithms ${Algorithms_SRCS} )

