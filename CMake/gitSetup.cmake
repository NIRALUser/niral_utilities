find_package(Git)
if(NOT GIT_FOUND)
  message(FATAL_ERROR "error: Install Git and try to re-configure")
endif()

#-----------------------------------------------------------------------------
# Git protocol option
#-----------------------------------------------------------------------------
# See https://github.com/Slicer/Slicer/blob/master/SuperBuild.cmake
option(${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
set(git_protocol "git")
if(NOT ${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL)
  set(git_protocol "http")
  # Verify that the global git config has been updated with the expected "insteadOf" option.
  function(_check_for_required_git_config_insteadof base insteadof)
  execute_process(
  COMMAND ${GIT_EXECUTABLE} config --global --get "url.${base}.insteadof"
  OUTPUT_VARIABLE output
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE error_code
  )
  if(error_code OR NOT "${output}" STREQUAL "${insteadof}")
    message(FATAL_ERROR
    "Since the ExternalProject modules doesn't provide a mechanism to customize the clone step by "
    "adding 'git config' statement between the 'git checkout' and the 'submodule init', it is required "
    "to manually update your global git config to successfully build ${CMAKE_PROJECT_NAME} with "
    "option ${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL set to FALSE. "
    "See http://na-mic.org/Mantis/view.php?id=2731"
    "\nYou could do so by running the command:\n"
    " ${GIT_EXECUTABLE} config --global url.\"${base}\".insteadOf \"${insteadof}\"\n")
    endif()
  endfunction()
  if("${ITK_VERSION_MAJOR}" LESS 4)
  _check_for_required_git_config_insteadof("http://itk.org/" "git://itk.org/")
  endif()
endif()

