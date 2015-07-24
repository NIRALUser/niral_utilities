project(NIRAL_UTILITIES)

cmake_minimum_required(VERSION 3.1)
cmake_policy(VERSION 3.1)

set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake
  ${CMAKE_CURRENT_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )
include(gitSetup)
include(niral_utilitiesExternalProject)
#-----------------------------------------------------------------------------

if(CMAKE_COMPILER_2005)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
endif()

if(${CMAKE_SOURCE_DIR} STREQUAL ${NIRAL_UTILITIES_SOURCE_DIR})
  set(LIBRARY_OUTPUT_PATH ${NIRAL_UTILITIES_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
  set(EXECUTABLE_OUTPUT_PATH ${NIRAL_UTILITIES_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")
endif()

set(ITK_IO_MODULES_USED
# IO Components
ITKCommon
ITKIOImageBase
ITKIOGDCM
ITKIOGIPL
ITKIOMeta
ITKIONIFTI
ITKIONRRD
ITKIOTIFF
ITKIOJPEG
ITKIOBMP
ITKIOPNG
ITKIOTransformInsightLegacy
ITKIOVTK
ITKIOSpatialObjects
ITKIOTransformBase
ITKIOHDF5
ITKIOTransformMatlab
ITKIOTransformHDF5
)

include_directories( Common )

option( BUILD_TESTING   "Build the testing tree" ON )
include(CTest)
include(ExternalData)

option(COMPILE_IMAGEMATH "Compile ImageMath." ON)
if(COMPILE_IMAGEMATH)
  subdirs(ImageMath)
endif()

option(COMPILE_IMAGESTAT "Compile ImageStat." OFF)
if(COMPILE_IMAGESTAT)
  subdirs(ImageStat)
endif()

option(COMPILE_CROPTOOLS "Compile CropTools." ON)
if(COMPILE_CROPTOOLS)
  subdirs(CropTools)
endif()

option(COMPILE_CONVERTITKFORMATS "Compile convertITKformats." OFF)
if(COMPILE_CONVERTITKFORMATS)
  subdirs(convertITKformats)
endif()

option(COMPILE_DWI_NIFTINRRDCONVERSION "Compile DWI_NiftiNrrdConversion." ON)
if(COMPILE_DWI_NIFTINRRDCONVERSION)
  subdirs(DWI_NiftiNrrdConversion)
endif()

option(COMPILE_POLYDATATRANSFORM "Compile PolydataTransform." ON)
if(COMPILE_POLYDATATRANSFORM)
  subdirs(PolydataTransform)
endif()

option(COMPILE_POLYDATAMERGE "Compile PolydataMerge." ON)
if(COMPILE_POLYDATAMERGE)
  subdirs(PolydataMerge)
endif()

option(COMPILE_CURVECOMPARE "Compile CurveCompare." OFF)
if(COMPILE_CURVECOMPARE)
  subdirs(CurveCompare)
endif()

option(COMPILE_MULTIATLASSEG "Compile MultiAtlasSeg." OFF)
if(COMPILE_MULTIATLASSEG)
  subdirs(MultiAtlasSeg)
endif()

option(COMPILE_MAURERDISTANCETRANSFORM "Compile MaurerDistanceTransform." OFF)
if(COMPILE_MAURERDISTANCETRANSFORM)
  subdirs(MaurerDistanceTransform)
endif()

option(COMPILE_POLYDATACOMPRESSION "Compile PolyDataCompression." OFF)
if(COMPILE_POLYDATACOMPRESSION)
  subdirs(PolyDataCompression)
endif()

niral_utilitiesExternalProject(
ITKTransformTools
GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/ITKTransformTools.git"
GIT_TAG "c7fb32d5d6f1cd6c0e57f1b5dfc0d1f58e371aae"
ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
DEPENDS ITK
INSTALL
)

niral_utilitiesExternalProject(
ReadImageHeader
GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/ReadImageHeader.git"
GIT_TAG "91484acae0e1fbb814ccccab329c77b1d1323362"
ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
DEPENDS SlicerExecutionModel ITK
INSTALL
)

niral_utilitiesExternalProject(
ImageNoise
GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/ImageNoise.git"
GIT_TAG "98e6a4056cc8d532d3e26201ba955342841fb515"
ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
DEPENDS SlicerExecutionModel ITK
INSTALL
)

niral_utilitiesExternalProject(
CombineDWIs
GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/CombineDWIs.git"
GIT_TAG "37d340b4b80a2b0735fbddcc450bf01d89b5b9bc"
ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
DEPENDS SlicerExecutionModel ITK
INSTALL
)

option(COMPILE_Unsupported "Compile Unsupported Tools." OFF)
if(COMPILE_Unsupported)
  subdirs(Unsupported)
endif()


