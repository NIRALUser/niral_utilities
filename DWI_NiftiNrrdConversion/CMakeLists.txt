cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

project(NiftiNrrdConversion)

find_package(GenerateCLP REQUIRED)
if(GenerateCLP_FOUND)
  include(${GenerateCLP_USE_FILE})
else()
   message(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
endif()

find_package(ITK COMPONENTS
  ITKCommon
  ITKIOImageBase
  ${ITK_IO_MODULES_USED}
)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
else()
  message(FATAL_ERROR "ITK not found. Please set ITK_DIR")
endif()

GENERATECLP(NiftiNrrdConversion.cxx NiftiNrrdConversion.xml)
add_executable(DWI_NiftiNrrdConversion NiftiNrrdConversion.cxx)
target_link_libraries(DWI_NiftiNrrdConversion ${ITK_LIBRARIES})

install( TARGETS DWI_NiftiNrrdConversion RUNTIME DESTINATION bin )
