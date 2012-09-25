
include(${AUTOSEG_SOURCE_DIR}/CMake/qtSetup.cmake)

INCLUDE_DIRECTORIES(
     	${CMAKE_CURRENT_SOURCE_DIR}
  	${CMAKE_CURRENT_BINARY_DIR}
	${QT_INCLUDE_DIR}
     )

LINK_LIBRARIES(
 ${QT_QT_LIBRARY}
)

SET(QtDisclaimer_SRCS 
  ../Disclaimer.cxx
  ../Disclaimer.h
  QtDisclaimer.cxx
  QtDisclaimer.h
)

SET(QtDisclaimer_GUI_SRCS 
   QtDisclaimerGUI.ui
)

IF(QT_WRAP_UI)
	QT_WRAP_UI(QtDisclaimer QtDisclaimer_HDRS  QtDisclaimer_SRCS  ${QtDisclaimer_GUI_SRCS} )
ENDIF(QT_WRAP_UI)

ADD_DEFINITIONS(-DQT_DLL)
ADD_DEFINITIONS(-DQT_THREAD_SUPPORT)


ADD_LIBRARY(QtDisclaimer ${QtDisclaimer_SRCS} )
