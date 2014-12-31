

find_package(ITK REQUIRED)
if(ITK_FOUND)
   include(${ITK_USE_FILE})
else(ITK_FOUND)
   message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
endif(ITK_FOUND)

include(${NIRAL_UTILITIES_SOURCE_DIR}/CMake/qtSetup.cmake)

include_directories(
     	#${NeuroLib_SOURCE_DIR}/Libraries/IO/Common
     	${NIRAL_UTILITIES_SOURCE_DIR}/Libraries/IO/XML
	${QT_INCLUDE_DIR}
     )

link_libraries(
  ITKIO
 ${QT_QT_LIBRARY}
)

add_definitions(-DQT_DLL)
add_definitions(-DQT_THREAD_SUPPORT)

set(IO_SRCS
  ${IO_SRCS}
  XML/XMLWriter.cxx
  XML/XMLReader.cxx
  XML/XMLIniIO.cxx
)

add_library(IO ${IO_SRCS} )

