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
ITKIOLSM
ITKIOStimulate
ITKIOBioRad
ITKIOMRC
ITKIOGE
MGHIO
)

include_directories( Common )

option( BUILD_TESTING   "Build the testing tree" ON )
if( BUILD_TESTING )
  include(CTest)
  include(ExternalData)
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

option(COMPILE_DTITensorNormalization "Compile DTITensorNormalization." OFF)
if(COMPILE_DTITensorNormalization)
  subdirs(DTITensorNormalization)
endif()

option(COMPILE_MAURERDISTANCETRANSFORM "Compile MaurerDistanceTransform." OFF)
if(COMPILE_MAURERDISTANCETRANSFORM)
  subdirs(MaurerDistanceTransform)
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

niral_utilitiesExternalProject(
NormalizeDWIGradientVectors
GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/NormalizeDWIGradientVectors.git"
GIT_TAG "e46bd5b7ca395ab7b290da621e329653287ac591"
ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
DEPENDS ITK
INSTALL
)

niral_utilitiesExternalProject(
3D4DImageConversion
GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/3D4DImageConversion.git"
GIT_TAG "5a53a98fce967a77c74d8e554d7231dc4c0ba695"
ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
DEPENDS SlicerExecutionModel ITK
INSTALL
)

niral_utilitiesExternalProject(
PolyDataCompression
GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/PolyDataCompression.git"
GIT_TAG "141b81712cb065fa890a1f49026495f2ea55fcc2"
ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
DEPENDS SlicerExecutionModel ITK VTK
INSTALL
)

option(COMPILE_Unsupported "Compile Unsupported Tools." OFF)
if(COMPILE_Unsupported)
  subdirs(Unsupported)
endif()

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH
  "Installation directory for header files")

if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/CMake/niral_utilities)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# The interesting stuff goes here
# ===============================
 
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE niral_utilities)
 
# Create the FooBarConfig.cmake and FooBarConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(niral_utilitiesConfig.cmake.in
  "${PROJECT_BINARY_DIR}/niral_utilitiesConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${NIRAL_UTILITIES_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(niral_utilitiesConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/niral_utilitiesConfig.cmake" @ONLY)
# ... for both
# configure_file(FooBarConfigVersion.cmake.in
#   "${PROJECT_BINARY_DIR}/FooBarConfigVersion.cmake" @ONLY)
 
# Install the FooBarConfig.cmake and FooBarConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/niral_utilitiesConfig.cmake"
  #"${PROJECT_BINARY_DIR}/FooBarConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
install(EXPORT niral_utilitiesTargets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)