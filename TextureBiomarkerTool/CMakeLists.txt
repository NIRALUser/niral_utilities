project(TextureBioMarkerTool)

cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

find_package(ITK)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
else()
  message(FATAL_ERROR
          "ITK not found. Please set ITK_DIR.")
endif()
list( FIND ITK_MODULES_ENABLED ITKReview FOUNDITKReview )
if( FOUNDITKReview EQUAL -1 )
  message(FATAL_ERROR "TextureBiomarkerTool: ITK needs to be compiled with Module_ITKReview set to ON")
endif()
add_executable(TextureBioMarkerTool TextureBioMarkerTool.cxx)
target_link_libraries(TextureBioMarkerTool ${ITK_LIBRARIES})

