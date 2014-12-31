# Locate Qt include paths and libraries

# This module defines
# QT_INCLUDE_DIR, where to find qt.h, etc.
# QT_LIBRARIES, the libraries to link against to use Qt.
# QT_DEFINITIONS, definitions to use when compiling code that uses Qt.
# QT_WRAP_CPP, If false, don't use QT_WRAP_CPP command.
# QT_WRAP_UI, If false, don't use QT_WRAP_UI command.
# QT_FOUND, If false, don't try to use Qt.

# also defined, but not for general use are
# QT_MOC_EXECUTABLE, where to find the moc tool.
# QT_UIC_EXECUTABLE, where to find the uic tool.
# QT_QT_LIBRARY, where to find the Qt library.
# QT_QTMAIN_LIBRARY, where to find the qtmain library. This is only required by Qt3 on Windows.


find_path(QT_INCLUDE_DIR qt.h
  $ENV{QTDIR}/include
  /usr/local/qt/include
  /usr/local/include
  /usr/include/qt3
  /usr/include/qt
  /usr/include
  C:/Progra~1/qt/include
  )

find_library(QT_QT_LIBRARY
  NAMES qt qt-mt qt-mt230nc
  PATHS
  $ENV{QTDIR}/lib
  /usr/local/qt/lib
  /usr/local/lib
  /usr/lib
  /usr/share/qt3/lib
  C:/Progra~1/qt/lib
  )

find_program(QT_MOC_EXECUTABLE moc
  $ENV{QTDIR}/bin C:/Progra~1/qt/bin
  )

find_program(QT_UIC_EXECUTABLE uic
  $ENV{QTDIR}/bin C:/Progra~1/qt/bin
  )


if(WIN32)
  find_library(QT_QTMAIN_LIBRARY qtmain
    $ENV{QTDIR}/lib C:/Progra~1/qt/lib
    DOC "This Library is only needed by and included with Qt3 on MSWindows. It should be NOTFOUND, undefined or IGNORE otherwise."
    )
endif(WIN32)


if(QT_MOC_EXECUTABLE)
  set( QT_WRAP_CPP "YES")
endif(QT_MOC_EXECUTABLE)

if(QT_UIC_EXECUTABLE)
  set( QT_WRAP_UI "YES")
endif(QT_UIC_EXECUTABLE)


if(QT_INCLUDE_DIR)
  if(QT_QT_LIBRARY)
    set( QT_LIBRARIES ${QT_LIBRARIES} ${QT_QT_LIBRARY} )
    set( QT_FOUND "YES" )
    set( QT_DEFINITIONS "")

    if(WIN32)
      if(QT_QTMAIN_LIBRARY)
        # for version 3
        set(QT_DEFINITIONS -DQT_DLL)
        set(QT_DEFINITIONS "-DQT_DLL -DQT_THREAD_SUPPORT -DNO_DEBUG")
        set(QT_LIBRARIES imm32.lib  ${QT_QT_LIBRARY} ${QT_QTMAIN_LIBRARY} )
        set(QT_LIBRARIES ${QT_LIBRARIES} winmm wsock32)
      else(QT_QTMAIN_LIBRARY)
        # for version 2
        set(QT_LIBRARIES imm32.lib ws2_32.lib  ${QT_QT_LIBRARY} )
      endif(QT_QTMAIN_LIBRARY)
    else(WIN32)
      set(QT_LIBRARIES ${QT_QT_LIBRARY} )
    endif(WIN32)

    # Backwards compatibility for CMake1.4 and 1.2
    set(QT_MOC_EXE ${QT_MOC_EXECUTABLE} )
    set(QT_UIC_EXE ${QT_UIC_EXECUTABLE} )

    if(UNIX)
      include( ${CMAKE_ROOT}/Modules/FindX11.cmake )
      if(X11_FOUND)
        set(QT_LIBRARIES ${QT_LIBRARIES} ${X11_LIBRARIES})
      endif(X11_FOUND)
      if(CMAKE_DL_LIBS)
        set(QT_LIBRARIES ${QT_LIBRARIES} ${CMAKE_DL_LIBS})
      endif(CMAKE_DL_LIBS)
    endif(UNIX)
    if(QT_QT_LIBRARY MATCHES "qt-mt")
      include( ${CMAKE_ROOT}/Modules/FindThreads.cmake )
      set(QT_LIBRARIES ${QT_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif(QT_QT_LIBRARY MATCHES "qt-mt")

  endif(QT_QT_LIBRARY)
else(QT_INCLUDE_DIR)
   message(FATAL_ERROR "QT library not found!\n" "Please go to http://www.ia.unc.edu/dev/tutorials/InstallLib")
endif(QT_INCLUDE_DIR)


mark_as_advanced(
  QT_INCLUDE_DIR
  QT_QT_LIBRARY
  QT_QTMAIN_LIBRARY
  QT_UIC_EXECUTABLE
  QT_MOC_EXECUTABLE
  QT_WRAP_CPP
  QT_WRAP_UI
  )

