
cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

project(TransformDeformationField)

find_package(GenerateCLP REQUIRED)
if(GenerateCLP_FOUND)
  include(${GenerateCLP_USE_FILE})
else(GenerateCLP_FOUND)
  message(FATAL_ERROR "GenerateCLP not found. Please set GenerateCLP_DIR.")
endif(GenerateCLP_FOUND)

find_package(ITK REQUIRED)
if(ITK_FOUND)
   include(${ITK_USE_FILE})
else(ITK_FOUND)
   message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
endif(ITK_FOUND)

if( ${ITK_VERSION_MAJOR} VERSION_GREATER 3 )
  message( FATAL_ERROR "TransformDeformationField needs ITKv3 to compile" )
endif()

set (CLP TransformDeformationField)

set ( ${CLP}_SOURCE ${CLP}.cpp )
generateclp(${CLP}_SOURCE ${CLP}.xml )

Set( dtiProcessFiles
dtiprocessFiles/deformationfieldio.h
dtiprocessFiles/deformationfieldio.cxx
dtiprocessFiles/itkHFieldToDeformationFieldImageFilter.h
dtiprocessFiles/itkHFieldToDeformationFieldImageFilter.txx
)

add_executable(${CLP} ${${CLP}_SOURCE} ${dtiProcessFiles} )

target_link_libraries (${CLP}
  ModuleDescriptionParser
  ITKIO
  ITKAlgorithms
  ITKBasicFilters
)


install( TARGETS ${CLP} RUNTIME DESTINATION bin )
