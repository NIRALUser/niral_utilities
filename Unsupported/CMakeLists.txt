message(WARNING "You are trying to build unsupported software. Beware, there might be errors")
option(COMPILE_IntensityRescaler "Compile IntensityRescaler" OFF)
if(COMPILE_IntensityRescaler)
  subdirs(IntensityRescaler)
endif()

