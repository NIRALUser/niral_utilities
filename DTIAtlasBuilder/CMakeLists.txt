PROJECT(DTIAtlasBuilder)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(VERSION 2.8)

FIND_PACKAGE(ITK REQUIRED)
IF (NOT ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF()

FIND_PACKAGE(GenerateCLP REQUIRED)
IF(NOT GenerateCLP_FOUND)
  MESSAGE(FATAL_ERROR "GenerateCLP not found. Please set GenerateCLP_DIR.")
ENDIF()


FIND_PACKAGE(Qt4 REQUIRED)
IF(NOT QT_USE_FILE)
  MESSAGE(FATAL_ERROR, "QT not found. Please set QT_DIR.")
ENDIF()


option(USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
set(git_protocol "git")
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "https")
endif()

include(ExternalProject)

MESSAGE( WARNING "Please verify that the following tools are installed on your system: ImageMath, ResampleDTIlogEuclidean, CropDTI, dtiprocess, BRAINSFit, AtlasWerks, dtiaverage, DTI-Reg, unu, MriWatcher")

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

