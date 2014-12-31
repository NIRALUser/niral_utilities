
option(COMPILE_DISPLAY "Compile Display library." ON)
mark_as_advanced(COMPILE_DISPLAY)

# if(USE_QT_LIB)
#   if(USE_ITK_LIB)
# 	option(COMPILE_IO "Compile IO library." ON)
# 	set(NeuroLib_COMPILE_IO ${COMPILE_IO})
# 	mark_as_advanced(COMPILE_IO)
#   endif(USE_ITK_LIB)
# endif(USE_QT_LIB)

#-----------------------------------------------------------------------------
# Compile disclaimer library
# subdirs(Disclaimer)
# subdirs(Display/Plot)
# subdirs(Datatypes)
# subdirs(IO/Qt)
# subdirs(IO/Common)
# subdirs(Algorithms)
# subdirs(Display/Common)
# subdirs(Display/Window2D)

#subdirs(Datatypes)
subdirs(IO)
subdirs(Algorithms)
subdirs(Display)

