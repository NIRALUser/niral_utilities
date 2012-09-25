
INCLUDE_DIRECTORIES(
 ${FLTK_INCLUDE_DIR}
 ${CMAKE_CURRENT_BINARY_DIR}
 ${CMAKE_CURRENT_SOURCE_DIR}
)


SET(FltkDisclaimer_UI
  FltkDisclaimerGUI.fl
)



SET(FltkDisclaimer_SRCS
  ../Disclaimer.h
  ../Disclaimer.cxx
  FltkDisclaimer.cxx
  FltkDisclaimer.h
  FltkDisclaimerGUIBase.h
)


FLTK_WRAP_UI(FltkDisclaimer ${FltkDisclaimer_UI})
ADD_LIBRARY(FltkDisclaimer ${FltkDisclaimer_SRCS} ${FltkDisclaimer_FLTK_UI_SRCS})
