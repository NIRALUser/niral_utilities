project(convertITKformats)

cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

find_package(ITK REQUIRED)
if(ITK_FOUND)
   include(${ITK_USE_FILE})
else(ITK_FOUND)
   message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
endif(ITK_FOUND)

if( ${ITK_VERSION_MAJOR} VERSION_GREATER 3)
  if( NOT ITKV3_COMPATIBILITY )
    message(FATAL_ERROR "convertITKformats: ITK needs to be compiled with ITKV3_COMPATIBILITY set to ON")
  endif()
endif()

add_executable( convertITKformats convertITKformats.cxx)
add_executable( convertITKformats4D convertITKformats4D.cxx)
add_executable( convertTextToImage convertTextToImage.cxx)

target_link_libraries( convertITKformats ${ITK_LIBRARIES} )
target_link_libraries( convertITKformats4D ${ITK_LIBRARIES} )
target_link_libraries( convertTextToImage ${ITK_LIBRARIES} )

install( TARGETS convertITKformats convertITKformats4D convertTextToImage RUNTIME DESTINATION bin )
