project(MultiAtlasSeg)

cmake_minimum_required(VERSION 3.5)

find_package(ITK REQUIRED)
if(ITK_FOUND)
   include(${ITK_USE_FILE})
else()
   message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
endif()

list( FIND ITK_MODULES_ENABLED ITKReview position)
if( "${position}" EQUAL "-1" )
  message( FATAL_ERROR "MultiAtlasSeg needs ITKv4 compiled with ITKReview enabled" )
endif()

add_executable(MultiAtlasSeg MultiAtlasSeg.cxx)
target_link_libraries(MultiAtlasSeg ${ITK_LIBRARIES})

install( TARGETS MultiAtlasSeg RUNTIME DESTINATION bin )
