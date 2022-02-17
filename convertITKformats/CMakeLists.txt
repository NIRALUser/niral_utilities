project(convertITKformats)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(ITK 5.3 REQUIRED COMPONENTS
  ITKCommon
  ITKDICOMParser
  ITKIOImageBase
  ITKImageFilterBase
  ITKImageIntensity
  ${ITK_IO_MODULES_USED}
  )
include(${ITK_USE_FILE})

#
# To remove the requirement on ITK built with ITKV3_COMPATIBILITY option enabled.
# This project includes its own copy of itkDICOMSeriesFileNames.h and
# itkDICOMSeriesFileNames.cxx:
#
# * to avoid conflicts when building against ITK build having ITKV3_COMPATIBILITY
#   set to ON. The fileshave been renamed to itkNiralDICOMSeriesFileNames and
#   namespace updated from "itk::" to "itk::niral::".
#
# * to find to which revision of ITK the files correspond, lookup the commit messages
#   associated with itkNiralDICOMSeriesFileNames.(h|cxx) from this project history.
#
# * instead of expecting ITKDeprecated module, only ITKDICOMParser module is required.
#   Since this module is enabled by default (it belongs on ITK_BUILD_DEFAULT_MODULES),
#   and is not deprecated this is a reasonable expectation.
#

add_definitions(-DDICOM_SUPPORT="True")

niral_add_executable(
  NAME convertITKformats
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  ADDITIONAL_SRCS itkNiralDICOMSeriesFileNames.cxx
  )
niral_add_executable(
  NAME convertITKformats4D
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  ADDITIONAL_SRCS itkNiralDICOMSeriesFileNames.cxx
  )
niral_add_executable(
  NAME convertTextToImage
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  )


install( TARGETS convertITKformats 
  EXPORT niral_utilitiesTargets
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
  )

set_property(GLOBAL APPEND PROPERTY niral_utilities_LIBRARIES convertITKformats)
