cmake_minimum_required(VERSION 2.6)

file(GLOB LIBRARY_SOURCE *.cxx *.txx *.h)

include_directories(${PROJECT_SOURCE_DIR}/Library)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(library ${LIBRARY_SOURCE})