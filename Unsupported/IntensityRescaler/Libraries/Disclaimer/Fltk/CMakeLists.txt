
include_directories(
 ${FLTK_INCLUDE_DIR}
 ${CMAKE_CURRENT_BINARY_DIR}
 ${CMAKE_CURRENT_SOURCE_DIR}
)


set(FltkDisclaimer_UI
  FltkDisclaimerGUI.fl
)



set(FltkDisclaimer_SRCS
  ../Disclaimer.h
  ../Disclaimer.cxx
  FltkDisclaimer.cxx
  FltkDisclaimer.h
  FltkDisclaimerGUIBase.h
)


fltk_wrap_ui(FltkDisclaimer ${FltkDisclaimer_UI})
add_library(FltkDisclaimer ${FltkDisclaimer_SRCS} ${FltkDisclaimer_FLTK_UI_SRCS})
