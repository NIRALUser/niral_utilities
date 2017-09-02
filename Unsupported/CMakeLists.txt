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

