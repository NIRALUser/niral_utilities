project(NIRAL_UTILITIES)

set(NIRAL_UTILITIES_MAJOR_VERSION 2)
set(NIRAL_UTILITIES_MINOR_VERSION 0)
set(NIRAL_UTILITIES_PATCH_VERSION 7)
set(NIRAL_UTILITIES_VERSION
  ${NIRAL_UTILITIES_MAJOR_VERSION}.${NIRAL_UTILITIES_MINOR_VERSION}.${NIRAL_UTILITIES_PATCH_VERSION})

cmake_minimum_required(VERSION 3.5)

#-----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake
  ${CMAKE_CURRENT_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )
include(gitSetup)
include(niral_add_executable)
include(niral_utilitiesExternalProject)

#-----------------------------------------------------------------------------
if(CMAKE_COMPILER_2005)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
endif()

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
#-----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Debug' as none was specified.")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  mark_as_advanced(CMAKE_BUILD_TYPE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

#-----------------------------------------------------------------------------
if(${CMAKE_SOURCE_DIR} STREQUAL ${NIRAL_UTILITIES_SOURCE_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${NIRAL_UTILITIES_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${NIRAL_UTILITIES_BINARY_DIR}/lib CACHE PATH "Single output directory for building all shared libraries.")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${NIRAL_UTILITIES_BINARY_DIR}/lib CACHE PATH "Single output directory for building all static libraries.")
endif()

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

#-----------------------------------------------------------------------------
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

if(NOT VTK_RENDERING_BACKEND)
  set(VTK_RENDERING_BACKEND OpenGL)
endif()

# Describe property
define_property(GLOBAL PROPERTY niral_utilities_LIBRARIES
    BRIEF_DOCS "Defined libraries"
    FULL_DOCS "Defined libraries")
# Initialize property
set_property(GLOBAL PROPERTY niral_utilities_LIBRARIES "")

include_directories( Common )

option( BUILD_TESTING   "Build the testing tree" ON )
if( BUILD_TESTING )
  include(CTest)
  include(ExternalData)
endif()

find_package(ITK)

if(ITK_VERSION)
  string(REPLACE "." ";" VERSION_LIST ${ITK_VERSION})
  list(GET VERSION_LIST 0 ITK_VERSION_MAJOR)
  list(GET VERSION_LIST 1 ITK_VERSION_MINOR)
  list(GET VERSION_LIST 2 ITK_VERSION_PATCH)

  if(ITK_VERSION_MAJOR GREATER 4 OR ITK_VERSION_MAJOR EQUAL 4 AND ITK_VERSION_MINOR GREATER 12)
    message(STATUS "USING compilation flag -std=c++11")
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FAGS} -std=c++11)
  endif()
endif()

option(COMPILE_IMAGEMATH "Compile ImageMath." ON)
if(COMPILE_IMAGEMATH)
  add_subdirectory(ImageMath)
endif()

option(COMPILE_IMAGESTAT "Compile ImageStat." OFF)
if(COMPILE_IMAGESTAT)
  add_subdirectory(ImageStat)
endif()

option(COMPILE_CROPTOOLS "Compile CropTools." ON)
if(COMPILE_CROPTOOLS)
  add_subdirectory(CropTools)
endif()

option(COMPILE_CONVERTITKFORMATS "Compile convertITKformats." OFF)
if(COMPILE_CONVERTITKFORMATS)
  add_subdirectory(convertITKformats)
endif()

option(COMPILE_DWI_NIFTINRRDCONVERSION "Compile DWI_NiftiNrrdConversion." ON)
if(COMPILE_DWI_NIFTINRRDCONVERSION)
  add_subdirectory(DWI_NiftiNrrdConversion)
endif()

option(COMPILE_POLYDATATRANSFORM "Compile PolydataTransform." ON)
if(COMPILE_POLYDATATRANSFORM)
  add_subdirectory(PolydataTransform)
endif()

option(COMPILE_POLYDATAMERGE "Compile PolydataMerge." ON)
if(COMPILE_POLYDATAMERGE)
  add_subdirectory(PolydataMerge)
endif()

option(COMPILE_MULTIATLASSEG "Compile MultiAtlasSeg." OFF)
if(COMPILE_MULTIATLASSEG)
  add_subdirectory(MultiAtlasSeg)
endif()

option(COMPILE_DTITensorNormalization "Compile DTITensorNormalization." OFF)
if(COMPILE_DTITensorNormalization)
  add_subdirectory(DTITensorNormalization)
endif()

option(COMPILE_MAURERDISTANCETRANSFORM "Compile MaurerDistanceTransform." OFF)
if(COMPILE_MAURERDISTANCETRANSFORM)
  add_subdirectory(MaurerDistanceTransform)
endif()

niral_utilitiesExternalProject(
  ITKTransformTools
  GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/ITKTransformTools.git"
  GIT_TAG "c7fb32d5d6f1cd6c0e57f1b5dfc0d1f58e371aae"
  ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
  DEPENDS ITK
  )

niral_utilitiesExternalProject(
  ReadImageHeader
  GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/ReadImageHeader.git"
  GIT_TAG "91484acae0e1fbb814ccccab329c77b1d1323362"
  ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
  DEPENDS SlicerExecutionModel ITK
  )

niral_utilitiesExternalProject(
  ImageNoise
  GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/ImageNoise.git"
  GIT_TAG "687531cc4ec659f93465a013ef2da153513e16db"
  ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
  DEPENDS SlicerExecutionModel ITK
  )

niral_utilitiesExternalProject(
  CombineDWIs
  GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/CombineDWIs.git"
  GIT_TAG "37d340b4b80a2b0735fbddcc450bf01d89b5b9bc"
  ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
  DEPENDS SlicerExecutionModel ITK
  )

niral_utilitiesExternalProject(
  NormalizeDWIGradientVectors
  GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/NormalizeDWIGradientVectors.git"
  GIT_TAG "e46bd5b7ca395ab7b290da621e329653287ac591"
  ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
  DEPENDS ITK
  )

niral_utilitiesExternalProject(
  3D4DImageConversion
  GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/3D4DImageConversion.git"
  GIT_TAG "5a53a98fce967a77c74d8e554d7231dc4c0ba695"
  ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
  DEPENDS SlicerExecutionModel ITK
  )

niral_utilitiesExternalProject(
  PolyDataCompression
  GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/PolyDataCompression.git"
  GIT_TAG 4650c5470d5dd621a4e421b1611c415d1c1ffb9b
  ADDITIONAL_OPTIONS
    -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
    -DVTK_RENDERING_BACKEND:STRING=${VTK_RENDERING_BACKEND}
    -DBUILD_TESTING:BOOL=${BUILD_TESTING}
  DEPENDS SlicerExecutionModel ITK VTK
  )

# CurveCompare: It compares two fibers file with different metrics: euclidean
# distance (closest points), tangent metric, curvature metric, overlap, probabilistic
# overlap. Works only with ITK4.
niral_utilitiesExternalProject(
  CurveCompare
  GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/CurveCompare.git"
  GIT_TAG 52eaf02611b12f771f1c82c93bbb760945946de4
  ADDITIONAL_OPTIONS -DCMAKE_INSTALL_PREFIX:PATH=${NIRAL_UTILITIES_BINARY_DIR}
  DEPENDS SlicerExecutionModel ITK VTK
  )

option(COMPILE_Unsupported "Compile Unsupported Tools." OFF)
if(COMPILE_Unsupported)
  add_subdirectory(Unsupported)
endif()

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






### This is the config file for niral_utilities
# Create the FooBarConfig.cmake and FooBarConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
get_property(CONF_LIBRARIES GLOBAL PROPERTY niral_utilities_LIBRARIES)

configure_file(niral_utilitiesConfig.cmake.in
  "${PROJECT_BINARY_DIR}/niral_utilitiesConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${NIRAL_UTILITIES_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(niral_utilitiesConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/niral_utilitiesConfig.cmake" @ONLY)

install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/niral_utilitiesConfig.cmake"
  #"${PROJECT_BINARY_DIR}/FooBarConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
if(COMPILE_POLYDATAMERGE OR COMPILE_POLYDATATRANSFORM)
 install(EXPORT niral_utilitiesTargets DESTINATION
"${INSTALL_CMAKE_DIR}" COMPONENT dev)
endif()

export(TARGETS ${CONF_LIBRARIES}
  FILE "${PROJECT_BINARY_DIR}/niral_utilitiesTargets.cmake"
  )

