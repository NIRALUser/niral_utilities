cmake_minimum_required(VERSION 2.6)

project(polydatamerge)

find_package(ITK COMPONENTS
  ITKCommon
  ITKIOImageBase
  ${ITK_IO_MODULES_USED}
)
include(${ITK_USE_FILE})

find_package(VTK COMPONENTS
    vtkCommonCore
    vtkFiltersCore
    vtkIOLegacy
)
include(${VTK_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR})


find_package(GenerateCLP REQUIRED)
if(GenerateCLP_FOUND)
  include(${GenerateCLP_USE_FILE})
else()
  message(FATAL_ERROR
          "GenerateCLP not found. Please set GenerateCLP_DIR.")
endif()

set(src polydatamerge.cpp)
GENERATECLP(src polydatamerge.xml)

add_executable(polydatamerge  ${src}  )
target_link_libraries(polydatamerge  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

install( TARGETS polydatamerge RUNTIME DESTINATION bin )
