

include_directories(
      ${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Datatypes
      ${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Algorithms
      ${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/IO/Common
      ${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Common
      ${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/Display/Plot
      ${ITK_SOURCE_DIR}/Utilities/zlib
     )

link_libraries(
 ${ITK_LIBRARIES}
)


#set(Display_MOC_SRCS
#)

#if(QT_WRAP_CPP)
#  qt_wrap_cpp(Display Displaymoced_SRCS ${Display_MOC_SRCS})
#endif()


set(Algorithms_SRCS
  ${Algorithms_SRCS}
  #${Algorithmsmoced_SRCS}
  ImageHistogram.cxx
  ImageIntensityNormalizer.cxx
  ImageCombinaison.cxx
)

if(COMPILE_IntensityRescaler_GUI)
  add_definitions(-DQT_DLL)
  add_definitions(-DQT_THREAD_SUPPORT)

  list( APPEND Algorithms_SRCS
  ImageRotation.cxx
  Rotation3D.cxx
  )
  link_libraries(
   ${QT_QT_LIBRARY}
   ${QT_LIBRARIES}
  )
  include_directories(
    ${QT_INCLUDE_DIR}
     )

endif()

add_library(Algorithms ${Algorithms_SRCS} )

