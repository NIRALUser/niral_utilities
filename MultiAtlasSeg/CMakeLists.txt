project(MultiAtlasSeg)

cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)


find_package(ITK REQUIRED)
if(ITK_FOUND)
   include(${ITK_USE_FILE})
else()
   message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
endif()

add_executable(MultiAtlasSeg MultiAtlasSeg.cxx)
target_link_libraries(MultiAtlasSeg ${ITK_LIBRARIES})

#install( TARGETS MultiAtlasSeg RUNTIME DESTINATION bin )
