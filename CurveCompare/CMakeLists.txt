PROJECT(CurveCompare)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(VERSION 2.8)

FIND_PACKAGE(ITK 4 REQUIRED)
IF (ITK_FOUND)
   INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR, "ITKv4 not found. Please set ITK_DIR.")
ENDIF (ITK_FOUND)

FIND_PACKAGE(GenerateCLP REQUIRED)
IF(GenerateCLP_FOUND)
  INCLUDE(${GenerateCLP_USE_FILE})
ELSE(GenerateCLP_FOUND)
  MESSAGE(FATAL_ERROR
          "GenerateCLP not found. Please set GenerateCLP_DIR.")
ENDIF(GenerateCLP_FOUND)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


IF( ${ITK_VERSION_MAJOR} VERSION_LESS 4 )
  MESSAGE( FATAL_ERROR "CurveCompare needs ITKv4 to compile" )
ENDIF()

include(ExternalProject)

option(USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
set(git_protocol "git")
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "http")
endif()

set(proj CurveCompare)
ExternalProject_Add(${proj}
    GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/CurveCompare.git"
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    CMAKE_ARGS 
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      -DGenerateCLP_DIR:PATH=${GenerateCLP_DIR}
      -DITK_DIR:PATH=${ITK_DIR}
      -DVTK_DIR:PATH=${VTK_DIR}
      -DEXECUTABLE_OUTPUT_PATH:PATH=${EXECUTABLE_OUTPUT_PATH}
    INSTALL_COMMAND ""
    DEPENDS ""
  )

