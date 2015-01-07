project(DTIAtlasBuilder)

cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

find_package(SlicerExecutionModel REQUIRED)
if(NOT SlicerExecutionModel_FOUND)
  message(FATAL_ERROR "SlicerExecutionModel not found. Please set
SlicerExecutionModel_DIR.")
endif()

find_package(ITK COMPONENTS
  ITKCommon
  ITKIOImageBase
  ${ITK_IO_MODULES_USED}
  REQUIRED
)

find_package(VTK COMPONENTS
          vtkCommonComputationalGeometry
          vtkCommonCore
          vtkCommonDataModel
          vtkCommonExecutionModel
          vtkCommonMath
          vtkCommonMisc
          vtkCommonSystem
          vtkCommonTransforms
          vtkFiltersCore
          vtkFiltersGeneral
          vtkFiltersModeling
          vtkFiltersSources
          vtkImagingCore
          vtkRenderingCore
          vtkRenderingOpenGL
          vtkTestingCore
          vtkTestingRendering
        )

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
      -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
      -DITK_DIR:PATH=${ITK_DIR}
      -DVTK_DIR:PATH=${VTK_DIR}
      -DQT_QMAKE_EXECUTABLE:PATH=${QT_QMAKE_EXECUTABLE}
      -DEXECUTABLE_OUTPUT_PATH:PATH=${EXECUTABLE_OUTPUT_PATH}
    INSTALL_COMMAND ""
    DEPENDS ""
  )

