project(CropTools)

if(CMAKE_COMPILER_2005)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
endif()

find_package(ITK 5.3 REQUIRED COMPONENTS
  ITKCommon
  ITKDiffusionTensorImage
  ITKIOImageBase
  ITKIOMeta
  ITKImageFilterBase
  ITKImageGrid
  ${ITK_IO_MODULES_USED}
  )
include(${ITK_USE_FILE})


foreach(appname IN ITEMS Crop4DImage CropDTI CropImage)
  niral_add_executable(
    NAME ${appname}
    TARGET_LIBRARIES ${ITK_LIBRARIES}
    )

  install( TARGETS ${appname} 
    EXPORT niral_utilitiesTargets
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
  )
# Add all targets to the build-tree export set

set_property(GLOBAL APPEND PROPERTY niral_utilities_LIBRARIES ${appname})

endforeach()
