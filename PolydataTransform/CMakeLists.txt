cmake_minimum_required(VERSION 2.6)

project(polydatatransform)

find_package(GenerateCLP REQUIRED)
if(GenerateCLP_FOUND)
  include(${GenerateCLP_USE_FILE})
else()
  message(FATAL_ERROR
          "GenerateCLP not found. Please set GenerateCLP_DIR.")
endif()

set(${PROJECT_NAME}_ITK_COMPONENTS
  ITKCommon
  ITKIOImageBase
  ITKIONRRD
  ITKImageFilterBase
  ITKIOTransformHDF5
  ITKIOTransformInsightLegacy
  ITKIOTransformMatlab
  ${ModuleDescriptionParser_ITK_COMPONENTS}
  )
find_package(ITK COMPONENTS ${${PROJECT_NAME}_ITK_COMPONENTS} REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/Library)
include_directories(${PROJECT_SOURCE_DIR})

add_subdirectory(Library)


set(src polydatatransform.cpp)
GENERATECLP(src polydatatransform.xml)

add_executable(polydatatransform  ${src} )
target_link_libraries( polydatatransform library ${ITK_LIBRARIES} ${VTK_LIBRARIES})

install( TARGETS polydatatransform RUNTIME DESTINATION bin )
