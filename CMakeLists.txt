project(NIRAL_UTILITIES)

cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake
  ${CMAKE_CURRENT_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )

include(CTest)

if(CMAKE_COMPILER_2005)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
endif(CMAKE_COMPILER_2005)

if(${CMAKE_SOURCE_DIR} STREQUAL ${NIRAL_UTILITIES_SOURCE_DIR})
  set(LIBRARY_OUTPUT_PATH ${NIRAL_UTILITIES_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
  set(EXECUTABLE_OUTPUT_PATH ${NIRAL_UTILITIES_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")
endif(${CMAKE_SOURCE_DIR} STREQUAL ${NIRAL_UTILITIES_SOURCE_DIR})

include_directories( Common )

option(COMPILE_CORREVAL "Compile CorrespondenceEvaluator." OFF)
if(COMPILE_CORREVAL)
  subdirs(CorrespondenceEvaluator)
endif(COMPILE_CORREVAL)

option(COMPILE_IMAGEMATH "Compile ImageMath." ON)
if(COMPILE_IMAGEMATH)
  subdirs(ImageMath)
endif(COMPILE_IMAGEMATH)

option(COMPILE_IMAGESTAT "Compile ImageStat." OFF)
if(COMPILE_IMAGESTAT)
  subdirs(ImageStat)
endif(COMPILE_IMAGESTAT)

option(COMPILE_CROPTOOLS "Compile CropTools." ON)
if(COMPILE_CROPTOOLS)
  subdirs(CropTools)
endif(COMPILE_CROPTOOLS)

option(COMPILE_CONVERTITKFORMATS "Compile convertITKformats." OFF)
if(COMPILE_CONVERTITKFORMATS)
  subdirs(convertITKformats)
endif(COMPILE_CONVERTITKFORMATS)

option(COMPILE_DWI_NIFTINRRDCONVERSION "Compile DWI_NiftiNrrdConversion." ON)
if(COMPILE_DWI_NIFTINRRDCONVERSION)
  subdirs(DWI_NiftiNrrdConversion)
endif(COMPILE_DWI_NIFTINRRDCONVERSION)

option(COMPILE_POLYDATATRANSFORM "Compile PolydataTransform." ON)
if(COMPILE_POLYDATATRANSFORM)
  subdirs(PolydataTransform)
endif(COMPILE_POLYDATATRANSFORM)

option(COMPILE_TRANSFORMDEFORMATIONFIELD "Compile TransformDeformationField." OFF)
if(COMPILE_TRANSFORMDEFORMATIONFIELD)
  subdirs(TransformDeformationField)
endif(COMPILE_TRANSFORMDEFORMATIONFIELD)

option(COMPILE_POLYDATAMERGE "Compile PolydataMerge." ON)
if(COMPILE_POLYDATAMERGE)
  subdirs(PolydataMerge)
endif(COMPILE_POLYDATAMERGE)

option(COMPILE_CURVECOMPARE "Compile CurveCompare." OFF)
if(COMPILE_CURVECOMPARE)
  subdirs(CurveCompare)
endif(COMPILE_CURVECOMPARE)

option(COMPILE_DTIAtlasBuilder "Compile DTIAtlasBuilder." OFF)
if(COMPILE_DTIAtlasBuilder)
  subdirs(DTIAtlasBuilder)
endif(COMPILE_DTIAtlasBuilder)

option(COMPILE_MULTIATLASSEG "Compile MultiAtlasSeg." OFF)
if(COMPILE_MULTIATLASSEG)
  subdirs(MultiAtlasSeg)
endif(COMPILE_MULTIATLASSEG)

option(COMPILE_TEXTUREBIOMARKERTOOL "Compile TextureBiomarkerTool." OFF)
if(COMPILE_TEXTUREBIOMARKERTOOL)
  subdirs(TextureBiomarkerTool)
endif(COMPILE_TEXTUREBIOMARKERTOOL)

option(COMPILE_DMDBIOMARKERTOOL "Compile DMDBiomarkerTool." OFF)
if(COMPILE_DMDBIOMARKERTOOL)
  subdirs(DMDBiomarkerTool)
endif(COMPILE_DMDBIOMARKERTOOL)




