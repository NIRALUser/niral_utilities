
include(CMakeParseArguments)
include(ExternalProject)


function( niral_utilitiesExternalProject projectname)
  set(options INSTALL)
  set(oneValueArgs GIT_REPOSITORY GIT_TAG)
  set(multiValueArgs ADDITIONAL_OPTIONS DEPENDS)
  cmake_parse_arguments(_ep "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(proj ${projectname})

  option(COMPILE_${proj} "Compile ${proj}" OFF)
  if(NOT COMPILE_${proj})
     return()
  endif()

  message(STATUS "Adding external project: ${projectname}")

  set(install_arg)
  if(NOT _ep_INSTALL)
    set(install_arg INSTALL_COMMAND ${CMAKE_COMMAND} -E echo)
  endif()

  foreach(var ${_ep_DEPENDS})
    find_package(${var} REQUIRED)
    include(${${var}_USE_FILE})
    list(APPEND _ep_ADDITIONAL_OPTIONS -D${var}_DIR:PATH=${${var}_DIR})
  endforeach()

  set( ${proj}_REPOSITORY ${_ep_GIT_REPOSITORY})
  set( ${proj}_GIT_TAG ${_ep_GIT_TAG} )

  if(APPLE)
    list(APPEND _ep_ADDITIONAL_OPTIONS
      -DCMAKE_BUNDLE_OUTPUT_DIRECTORY:PATH=${CMAKE_BUNDLE_OUTPUT_DIRECTORY}
      )
  endif()

  if(DEFINED CMAKE_BUILD_TYPE)
    list(APPEND _ep_ADDITIONAL_OPTIONS
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      )
  endif()

  # Disable the "You are in 'detached HEAD' state." warning.
  set(git_config_arg)
  if(CMAKE_VERSION VERSION_GREATER "3.7.2")
    set(git_config_arg GIT_CONFIG "advice.detachedHead=false")
  endif()

  ExternalProject_Add(${projectname}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    ${git_config_arg}
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DCMAKE_SKIP_RPATH:BOOL=${CMAKE_SKIP_RPATH}
      -DCMAKE_MODULE_PATH:PATH=${CMAKE_MODULE_PATH}
      -DCMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
      -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
      -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
      -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
      -DCMAKE_MODULE_LINKER_FLAGS:STRING=${CMAKE_MODULE_LINKER_FLAGS}
      -DCMAKE_GENERATOR:STRING=${CMAKE_GENERATOR}
      -DCMAKE_EXTRA_GENERATOR:STRING=${CMAKE_EXTRA_GENERATOR}
      -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
      -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
      -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
      ${_ep_ADDITIONAL_OPTIONS}
    ${install_arg}
  )
endfunction()
