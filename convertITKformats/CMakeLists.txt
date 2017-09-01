project(convertITKformats)

cmake_minimum_required(VERSION 3.5)

message(WARNING "convertITKformats requires the ITK module ITKDeprecated.")

find_package(ITK COMPONENTS
  ITKCommon
  ITKDeprecated
  ITKIOImageBase
  ITKImageFilterBase
  ITKImageIntensity
  ${ITK_IO_MODULES_USED}
)

if(ITK_FOUND)
   include(${ITK_USE_FILE})
else()
   message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
endif()

if( ${ITK_VERSION_MAJOR} VERSION_GREATER 3)
  if( NOT ITKV3_COMPATIBILITY )
    message(WARNING "convertITKformats: To support DICOM, ITK needs to be compiled with ITKV3_COMPATIBILITY set to ON")
  else()
    list( FIND ITK_MODULES_ENABLED ITKIODCMTK position )
    if( ${position} GREATER -1 ) #not found: ${position}==-1
      ADD_DEFINITIONS(-DDICOM_SUPPORT="True")
    else()
      message(WARNING "convertITKformats: To support DICOM, ITK needs to be compiled with ITKIODCMTK")
    endif()
  endif()
endif()

add_executable( convertITKformats convertITKformats.cxx)
add_executable( convertITKformats4D convertITKformats4D.cxx)
add_executable( convertTextToImage convertTextToImage.cxx)

target_link_libraries( convertITKformats ${ITK_LIBRARIES} )
target_link_libraries( convertITKformats4D ${ITK_LIBRARIES} )
target_link_libraries( convertTextToImage ${ITK_LIBRARIES} )

install( TARGETS convertITKformats convertITKformats4D convertTextToImage RUNTIME DESTINATION bin )
