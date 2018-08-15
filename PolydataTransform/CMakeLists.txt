cmake_minimum_required(VERSION 3.5)

#project(polydatatransform)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED COMPONENTS
  vtkCommonCore
  vtkFiltersCore
  vtkIOLegacy
  vtkIOXML)
include(${VTK_USE_FILE})

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Library)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(Library)

set(LOCAL_NAME polydatatransform)
set(LOCAL_LIBRARIES 
	${ITK_LIBRARIES}
	${VTK_LIBRARIES}
  library
  )

SEMMacroBuildCLI(
  NAME ${LOCAL_NAME}
  EXECUTABLE_ONLY
  TARGET_LIBRARIES ${LOCAL_LIBRARIES}
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
  INSTALL_RUNTIME_DESTINATION ${INSTALL_RUNTIME_DESTINATION}
  INSTALL_LIBRARY_DESTINATION ${INSTALL_LIBRARY_DESTINATION}
  INSTALL_ARCHIVE_DESTINATION ${INSTALL_ARCHIVE_DESTINATION}
  )

file(GLOB LIBRARY_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Library/*.txx ${CMAKE_CURRENT_SOURCE_DIR}/Library/*.hxx ${CMAKE_CURRENT_SOURCE_DIR}//Library/*.h)

set_target_properties(polydatatransform PROPERTIES
  PUBLIC_HEADER "${LIBRARY_HEADERS}"
  )

install( TARGETS polydatatransform 
  EXPORT niral_utilitiesTargets
  RUNTIME DESTINATION bin 
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
  )
# Add all targets to the build-tree export set

set_property(GLOBAL APPEND PROPERTY niral_utilities_LIBRARIES polydatatransform)