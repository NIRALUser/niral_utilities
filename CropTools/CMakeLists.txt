project(CropTools)

cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

if(CMAKE_COMPILER_2005)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
endif()

find_package(ITK REQUIRED)
if(ITK_FOUND)
   include(${ITK_USE_FILE})
else()
   message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
endif()


add_executable(Crop4DImage Crop4DImage.cxx )
target_link_libraries( Crop4DImage ${ITK_LIBRARIES} )
add_executable(CropDTI CropDTI.cxx )
target_link_libraries( CropDTI ${ITK_LIBRARIES} )
add_executable(CropImage CropImage.cxx )
target_link_libraries( CropImage ${ITK_LIBRARIES} )

install( TARGETS Crop4DImage CropDTI CropImage RUNTIME DESTINATION bin )
