cmake_minimum_required(VERSION 3.5)

project(DTITensorNormalization)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable( DTITensorNormalization DTITensorNormalization.cxx)
target_link_libraries(DTITensorNormalization ${ITK_LIBRARIES} )

install( TARGETS DTITensorNormalization RUNTIME DESTINATION bin )

