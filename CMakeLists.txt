#-----------------------------------------------------------------------------
# Define cmake build policies that affect compilation and linkage default behaviors
#
# Set the NEWEST_VALIDATED_POLICIES_VERSION string to the newest cmake version
# policies that provide successful builds. By setting NEWEST_VALIDATED_POLICIES_VERSION
# to a value greater than the oldest policies, all policies between
# OLDEST_VALIDATED_POLICIES_VERSION and CMAKE_VERSION (used for this build)
# are set to their NEW behaivor, thereby suppressing policy warnings related to policies
# between the OLDEST_VALIDATED_POLICIES_VERSION and CMAKE_VERSION.
#
# CMake versions greater than the NEWEST_VALIDATED_POLICIES_VERSION policies will
# continue to generate policy warnings "CMake Warning (dev)...Policy CMP0XXX is not set:"
#
set(OLDEST_VALIDATED_POLICIES_VERSION "3.16.3")
set(NEWEST_VALIDATED_POLICIES_VERSION "3.19.7")
cmake_minimum_required(VERSION ${OLDEST_VALIDATED_POLICIES_VERSION}...${NEWEST_VALIDATED_POLICIES_VERSION} FATAL_ERROR)

#-----------------------------------------------------------------------------
# Enable C++14
#-----------------------------------------------------------------------------
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message(STATUS "${_msg} - C++${CMAKE_CXX_STANDARD}")
if(NOT CMAKE_CXX_STANDARD MATCHES "^(14|17|20)$")
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be set to 14, 17, or 20")
endif()


set(NIRAL_UTILITIES_MAJOR_VERSION 2)
set(NIRAL_UTILITIES_MINOR_VERSION 0)
set(NIRAL_UTILITIES_PATCH_VERSION 7)
set(NIRAL_UTILITIES_VERSION
  ${NIRAL_UTILITIES_MAJOR_VERSION}.${NIRAL_UTILITIES_MINOR_VERSION}.${NIRAL_UTILITIES_PATCH_VERSION})

project(NIRAL_UTILITIES 
    VERSION ${NIRAL_UTILITIES_MAJOR_VERSION}.${NIRAL_UTILITIES_MINOR_VERSION}.${NIRAL_UTILITIES_PATCH_VERSION}
    LANGUAGES CXX C)

#-----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake
  ${CMAKE_CURRENT_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )
include(gitSetup)
include(niral_add_executable)
include(niral_utilitiesExternalProject)

if(NOT COMMAND SETIFEMPTY)
  macro(SETIFEMPTY)
    set(KEY ${ARGV0})
    set(VALUE ${ARGV1})
    if(NOT ${KEY})
      set(${ARGV})
    endif()
  endmacro()
endif()

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
  set(VTK_RENDERING_BACKEND OpenGL2)
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

find_package(ITK 5.3 REQUIRED)

setifempty(INSTALL_RUNTIME_DESTINATION bin)
setifempty(INSTALL_LIBRARY_DESTINATION lib)
setifempty(INSTALL_ARCHIVE_DESTINATION lib)

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
set(CONF_INCLUDE_DIRS "\${niral_utilities_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(niral_utilitiesConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/niral_utilitiesConfig.cmake" @ONLY)

install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/niral_utilitiesConfig.cmake"
  #"${PROJECT_BINARY_DIR}/FooBarConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree

install(EXPORT niral_utilitiesTargets DESTINATION
"${INSTALL_CMAKE_DIR}" COMPONENT dev)


export(TARGETS ${CONF_LIBRARIES}
  FILE "${PROJECT_BINARY_DIR}/niral_utilitiesTargets.cmake"
  )

