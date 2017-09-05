

include_directories(
  #${NeuroLib_SOURCE_DIR}/Libraries/IO/Common
  ${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/IO/XML
	${QT_INCLUDE_DIR}
  )

link_libraries(
  ITKIO
  ${QT_QT_LIBRARY}
  ${QT_LIBRARIES}
  )

add_definitions(-DQT_DLL)
add_definitions(-DQT_THREAD_SUPPORT)


set(IO_SRCS
  ${IO_SRCS}
  XML/XMLWriter.cxx
  XML/XMLReader.cxx
  XML/XMLIniIO.cxx
  )

if( COMPILE_IntensityRescaler_GUI )
  include_directories(
  	${QT_INCLUDE_DIR}
    )
  link_libraries(
    ${QT_QT_LIBRARY}
    ${QT_LIBRARIES}
    )
  add_definitions(-DQT_DLL)
  add_definitions(-DQT_THREAD_SUPPORT)
endif()

add_library(IO ${IO_SRCS} )

