

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

include(${NIRAL_UTILITIES_SOURCE_DIR}/CMake/qtSetup.cmake)

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
)


#set(Display_MOC_SRCS
#)

#if(QT_WRAP_CPP)
#  qt_wrap_cpp(Display Displaymoced_SRCS ${Display_MOC_SRCS})
#endif(QT_WRAP_CPP)

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

