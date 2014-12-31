cmake_minimum_required(VERSION 2.6)
 
PROJECT(polydatamerge)


find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
 
find_package(VTK COMPONENTS
    vtkCommonCore
    vtkFiltersCore
    vtkIOLegacy
)
include(${VTK_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR})


FIND_PACKAGE(GenerateCLP REQUIRED)
IF(GenerateCLP_FOUND)
  INCLUDE(${GenerateCLP_USE_FILE})
ELSE(GenerateCLP_FOUND)
  MESSAGE(FATAL_ERROR
          "GenerateCLP not found. Please set GenerateCLP_DIR.")
ENDIF(GenerateCLP_FOUND)

SET(src polydatamerge.cpp)
GENERATECLP(src polydatamerge.xml)

add_executable(polydatamerge  ${src}  )
target_link_libraries(polydatamerge  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

install( TARGETS polydatamerge RUNTIME DESTINATION bin )
