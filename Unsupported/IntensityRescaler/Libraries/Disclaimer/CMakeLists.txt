if(COMPILE_IntensityRescaler_GUI)
  subdirs(Qt)
endif()

find_package(FLTK QUIET)
if(FLTK_FOUND)
    subdirs(Fltk)
else()
  message(WARNING, "FLTK not found. FLTK disclaimer will not be built. please set FLTK_DIR.")
endif()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

set(Disclaimer_SRCS
  Disclaimer.cxx
  Disclaimer.h
  )

add_library(Disclaimer ${Disclaimer_SRCS} )
