project(DTIAtlasBuilder)

cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

find_package(ITK REQUIRED)
if(NOT ITK_FOUND)
  message(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
endif()

find_package(GenerateCLP REQUIRED)
if(NOT GenerateCLP_FOUND)
  message(FATAL_ERROR "GenerateCLP not found. Please set GenerateCLP_DIR.")
endif()


find_package(Qt4 REQUIRED)
if(NOT QT_USE_FILE)
  message(FATAL_ERROR, "QT not found. Please set QT_DIR.")
endif()


option(USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
set(git_protocol "git")
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "https")
endif()

include(ExternalProject)

message( WARNING "Please verify that the following tools are installed on your system: ImageMath, ResampleDTIlogEuclidean, CropDTI, dtiprocess, BRAINSFit, AtlasWerks, dtiaverage, DTI-Reg, unu, MriWatcher")

set(proj DTIAtlasBuilder)
ExternalProject_Add(${proj}
    GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/DTIAtlasBuilder.git"
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      -DGenerateCLP_DIR:PATH=${GenerateCLP_DIR}
      -DITK_DIR:PATH=${ITK_DIR}
      -DQT_QMAKE_EXECUTABLE:PATH=${QT_QMAKE_EXECUTABLE}
      -DEXECUTABLE_OUTPUT_PATH:PATH=${EXECUTABLE_OUTPUT_PATH}
    INSTALL_COMMAND ""
    DEPENDS ""
  )

