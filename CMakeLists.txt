project(NIRAL_UTILITIES)

cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake
  ${CMAKE_CURRENT_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )


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

option(COMPILE_CORREVAL "Compile CorrespondenceEvaluator." OFF)
if(COMPILE_CORREVAL)
  subdirs(CorrespondenceEvaluator)
endif()

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

option(COMPILE_TRANSFORMDEFORMATIONFIELD "Compile TransformDeformationField." OFF)
if(COMPILE_TRANSFORMDEFORMATIONFIELD)
  subdirs(TransformDeformationField)
endif()

option(COMPILE_POLYDATAMERGE "Compile PolydataMerge." ON)
if(COMPILE_POLYDATAMERGE)
  subdirs(PolydataMerge)
endif()

option(COMPILE_CURVECOMPARE "Compile CurveCompare." OFF)
if(COMPILE_CURVECOMPARE)
  subdirs(CurveCompare)
endif()

option(COMPILE_DTIAtlasBuilder "Compile DTIAtlasBuilder." OFF)
if(COMPILE_DTIAtlasBuilder)
  subdirs(DTIAtlasBuilder)
endif()

option(COMPILE_MULTIATLASSEG "Compile MultiAtlasSeg." OFF)
if(COMPILE_MULTIATLASSEG)
  subdirs(MultiAtlasSeg)
endif()

option(COMPILE_TEXTUREBIOMARKERTOOL "Compile TextureBiomarkerTool." OFF)
if(COMPILE_TEXTUREBIOMARKERTOOL)
  subdirs(TextureBiomarkerTool)
endif()

option(COMPILE_DMDBIOMARKERTOOL "Compile DMDBiomarkerTool." OFF)
if(COMPILE_DMDBIOMARKERTOOL)
  subdirs(DMDBiomarkerTool)
endif()

option(COMPILE_Unsupported "Compile Unsupported Tools." OFF)
if(COMPILE_Unsupported)
  subdirs(Unsupported)
endif()


