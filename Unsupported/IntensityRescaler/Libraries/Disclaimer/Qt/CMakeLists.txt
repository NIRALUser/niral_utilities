
include(${AUTOSEG_SOURCE_DIR}/CMake/qtSetup.cmake)

include_directories(
     	${CMAKE_CURRENT_SOURCE_DIR}
  	${CMAKE_CURRENT_BINARY_DIR}
	${QT_INCLUDE_DIR}
     )

link_libraries(
 ${QT_QT_LIBRARY}
)

set(QtDisclaimer_SRCS
  ../Disclaimer.cxx
  ../Disclaimer.h
  QtDisclaimer.cxx
  QtDisclaimer.h
)

set(QtDisclaimer_GUI_SRCS
   QtDisclaimerGUI.ui
)

if(QT_WRAP_UI)
	qt_wrap_ui(QtDisclaimer QtDisclaimer_HDRS  QtDisclaimer_SRCS  ${QtDisclaimer_GUI_SRCS} )
endif(QT_WRAP_UI)

add_definitions(-DQT_DLL)
add_definitions(-DQT_THREAD_SUPPORT)


add_library(QtDisclaimer ${QtDisclaimer_SRCS} )
