cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

project(NiftiNrrdConversion)

find_package(ITK)
if(ITK_FOUND)
	include(${ITK_USE_FILE})
else(ITK_FOUND)
	message(FATAL_ERROR "ITK not found. Please set ITK_DIR")
endif(ITK_FOUND)

find_package(GenerateCLP REQUIRED)
if(GenerateCLP_FOUND)
  include(${GenerateCLP_USE_FILE})
else(GenerateCLP_FOUND)
   message(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
endif(GenerateCLP_FOUND)

GENERATECLP(NiftiNrrdConversion.cxx NiftiNrrdConversion.xml)
add_executable(DWI_NiftiNrrdConversion NiftiNrrdConversion.cxx)

target_link_libraries(DWI_NiftiNrrdConversion ${ITK_LIBRARIES})

install( TARGETS DWI_NiftiNrrdConversion RUNTIME DESTINATION bin )
