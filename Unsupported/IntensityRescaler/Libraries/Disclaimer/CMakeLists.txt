IF(USE_QT_LIB)
  SUBDIRS(Qt)
ENDIF(USE_QT_LIB)

FIND_PACKAGE(FLTK QUIET)
IF(FLTK_FOUND)
    SUBDIRS(Fltk)
ELSE(FLTK_FOUND)
  MESSAGE(WARNING, "FLTK not found. FLTK disclaimer will not be built. please set FLTK_DIR.")
ENDIF(FLTK_FOUND)

INCLUDE_DIRECTORIES(
     	${CMAKE_CURRENT_SOURCE_DIR}
  	${CMAKE_CURRENT_BINARY_DIR}
     )

SET(Disclaimer_SRCS 
  Disclaimer.cxx
  Disclaimer.h
)


ADD_LIBRARY(Disclaimer ${Disclaimer_SRCS} )
