project(CurveCompare)

cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

find_package(ITK 4 REQUIRED)
if(ITK_FOUND)
   include(${ITK_USE_FILE})
else(ITK_FOUND)
   message(FATAL_ERROR, "ITKv4 not found. Please set ITK_DIR.")
endif(ITK_FOUND)

find_package(GenerateCLP REQUIRED)
if(GenerateCLP_FOUND)
  include(${GenerateCLP_USE_FILE})
else(GenerateCLP_FOUND)
  message(FATAL_ERROR
          "GenerateCLP not found. Please set GenerateCLP_DIR.")
endif(GenerateCLP_FOUND)

find_package(VTK COMPONENTS
    vtkCommonCore
    vtkCommonExecutionModel
    vtkFiltersSources
    vtkIOLegacy
    vtkIOXML
    vtkInteractionWidgets
    vtkViewsInfovis
)
include(${VTK_USE_FILE})

if( ${ITK_VERSION_MAJOR} VERSION_LESS 4 )
  message( FATAL_ERROR "CurveCompare needs ITKv4 to compile" )
endif()

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

