
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(VERSION 2.8)

project(TransformDeformationField)

FIND_PACKAGE(GenerateCLP REQUIRED)
IF(GenerateCLP_FOUND)
  INCLUDE(${GenerateCLP_USE_FILE})
ELSE(GenerateCLP_FOUND)
  MESSAGE(FATAL_ERROR "GenerateCLP not found. Please set GenerateCLP_DIR.")
ENDIF(GenerateCLP_FOUND)

FIND_PACKAGE(ITK REQUIRED)
IF (ITK_FOUND)
   INCLUDE(${USE_ITK_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
ENDIF (ITK_FOUND)

IF( ${ITK_VERSION_MAJOR} VERSION_GREATER 3 )
  MESSAGE( FATAL_ERROR "TransformDeformationField needs ITKv3 to compile" )
ENDIF()

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
