message(WARNING "You are trying to build unsupported software. Beware, there might be errors")
option(COMPILE_IntensityRescaler "Compile IntensityRescaler" OFF)
if(COMPILE_IntensityRescaler)
  add_subdirectory(IntensityRescaler)
endif()

option(COMPILE_TRANSFORMDEFORMATIONFIELD "Compile TransformDeformationField." OFF)
if(COMPILE_TRANSFORMDEFORMATIONFIELD)
  add_subdirectory(TransformDeformationField)
endif()

option(COMPILE_CORREVAL "Compile CorrespondenceEvaluator." OFF)
if(COMPILE_CORREVAL)
  add_subdirectory(CorrespondenceEvaluator)
endif()

option(COMPILE_DMDBIOMARKERTOOL "Compile DMDBiomarkerTool." OFF)
if(COMPILE_DMDBIOMARKERTOOL)
  add_subdirectory(DMDBiomarkerTool)
endif()

option(COMPILE_TEXTUREBIOMARKERTOOL "Compile TextureBiomarkerTool." OFF)
if(COMPILE_TEXTUREBIOMARKERTOOL)
  add_subdirectory(TextureBiomarkerTool)
endif()

option(COMPILE_DTITensorNormalization "Compile DTITensorNormalization." OFF)
if(COMPILE_DTITensorNormalization)
  add_subdirectory(DTITensorNormalization)
endif()

option(COMPILE_DWI_NIFTINRRDCONVERSION "Compile DWI_NiftiNrrdConversion." ON)
if(COMPILE_DWI_NIFTINRRDCONVERSION)
  add_subdirectory(DWI_NiftiNrrdConversion)
endif()

niral_utilitiesExternalProject(
  3D4DImageConversion
  GIT_REPOSITORY "${git_protocol}://github.com/NIRALUser/3D4DImageConversion.git"
  GIT_TAG "5a53a98fce967a77c74d8e554d7231dc4c0ba695"
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
