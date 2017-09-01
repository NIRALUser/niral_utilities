project(CurveCompare)

cmake_minimum_required(VERSION 3.5)

find_package(GenerateCLP REQUIRED)
if(GenerateCLP_FOUND)
  include(${GenerateCLP_USE_FILE})
else()
  message(FATAL_ERROR
          "GenerateCLP not found. Please set GenerateCLP_DIR.")
endif()

find_package(ITK COMPONENTS
  ITKCommon
  ITKIOImageBase
  ITKImageAdaptors
  ITKImageFunction
  ITKImageStatistics
  ITKRegistrationCommon
  ITKSpatialObjects
  ITKTransform
  ${ITK_IO_MODULES_USED}
  REQUIRED
)
include(${ITK_USE_FILE})

find_package(VTK COMPONENTS
  vtkCommonComputationalGeometry
  vtkCommonCore
  vtkCommonDataModel
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

