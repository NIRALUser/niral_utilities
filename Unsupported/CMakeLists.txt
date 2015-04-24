message(WARNING "You are trying to build unsupported software. Beware, there might be errors")
option(COMPILE_IntensityRescaler "Compile IntensityRescaler" OFF)
if(COMPILE_IntensityRescaler)
  subdirs(IntensityRescaler)
endif()

option(COMPILE_CONVERTITKFORMATS "Compile convertITKformats." OFF)
if(COMPILE_CONVERTITKFORMATS)
  subdirs(convertITKformats)
endif()

option(COMPILE_TRANSFORMDEFORMATIONFIELD "Compile TransformDeformationField." OFF)
if(COMPILE_TRANSFORMDEFORMATIONFIELD)
  subdirs(TransformDeformationField)
endif()

option(COMPILE_CORREVAL "Compile CorrespondenceEvaluator." OFF)
if(COMPILE_CORREVAL)
  subdirs(CorrespondenceEvaluator)
endif()

option(COMPILE_DMDBIOMARKERTOOL "Compile DMDBiomarkerTool." OFF)
if(COMPILE_DMDBIOMARKERTOOL)
  subdirs(DMDBiomarkerTool)
endif()

option(COMPILE_TEXTUREBIOMARKERTOOL "Compile TextureBiomarkerTool." OFF)
if(COMPILE_TEXTUREBIOMARKERTOOL)
  subdirs(TextureBiomarkerTool)
endif()

