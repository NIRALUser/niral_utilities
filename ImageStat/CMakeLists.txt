cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

project(ImageStat)

find_package(ITK REQUIRED)
if(ITK_FOUND)
   include(${ITK_USE_FILE})
else()
   message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
endif()

add_executable( ImageStat ImageStat.cc)
target_link_libraries( ImageStat ${ITK_LIBRARIES} )

if( NOT ${ITK_USE_REVIEW_STATISTICS} )
  message(FATAL_ERROR "ITK has been compiled with ITK_USE_REVIEW_STATISTICS
set to OFF. ImageStat will fail to compile if this option is OFF. Please
recompile ITK with this option set to ON to compile ImageStat")
endif()

install( TARGETS ImageStat RUNTIME DESTINATION bin )
