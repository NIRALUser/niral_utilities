project(CropTools)

cmake_minimum_required(VERSION 3.5)

if(CMAKE_COMPILER_2005)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
endif()

find_package(ITK REQUIRED COMPONENTS
  ITKCommon
  ITKDiffusionTensorImage
  ITKIOImageBase
  ITKIOMeta
  ITKImageFilterBase
  ITKImageGrid
  ${ITK_IO_MODULES_USED}
  )
include(${ITK_USE_FILE})

add_executable(Crop4DImage Crop4DImage.cxx )
target_link_libraries( Crop4DImage ${ITK_LIBRARIES} )
add_executable(CropDTI CropDTI.cxx )
target_link_libraries( CropDTI ${ITK_LIBRARIES} )
add_executable(CropImage CropImage.cxx )
target_link_libraries( CropImage ${ITK_LIBRARIES} )

install( TARGETS Crop4DImage CropDTI CropImage RUNTIME DESTINATION bin )
