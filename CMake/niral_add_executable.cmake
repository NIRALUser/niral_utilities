
#
# This function builds a standalone NIRAL utility executable.
#
# Usage:
#
#   niral_add_executable(
#     NAME name
#     [SRC aname.cpp]
#     [ADDITIONAL_SRCS other1.cpp [other2.cpp [...]]]
#     [TARGET_LIBRARIES lib1 [lib2 [...]]]
#     )
#
# If only the NAME is specified, the function will lookup for source files like name.cpp, name.cxx,
# name.cc, ... (complete list of extensions considered are the one associated with add_executable CMake command)
#
# If the main executable source file has different name, the SRC allows to explicitly specify it.
#
# ADDITIONAL_SRCS allows to list any other source files that should be compiled in.
#
# TARGET_LIBRARIES allows to list all libraries to link the executable with.
#

function(niral_add_executable)
  set(options)
  set(oneValueArgs NAME SRC)
  set(multiValueArgs ADDITIONAL_SRCS TARGET_LIBRARIES)
  cmake_parse_arguments(MY "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(MY_UNPARSED_ARGUMENTS)
    message(AUTHOR_WARNING "Unparsed arguments given [${MY_UNPARSED_ARGUMENTS}]")
  endif()

  message(STATUS "Configuring executable: ${MY_NAME}")
  set(_src ${MY_NAME}.cxx)
  if(DEFINED MY_SRC)
    set(_src ${MY_SRC})
  endif()
  add_executable( ${MY_NAME} ${_src} ${MY_ADDITIONAL_SRCS})
  if(MY_TARGET_LIBRARIES)
    target_link_libraries(${MY_NAME} ${MY_TARGET_LIBRARIES})
  endif()
  install(TARGETS ${MY_NAME} RUNTIME DESTINATION bin)
endfunction()
