# project ImageMath
# Add standard remote object stores to user's
# configuration.
list(APPEND ExternalData_URL_TEMPLATES
  "http://slicer.kitware.com/midas3/api/rest?method=midas.bitstream.download&checksum=%(hash)"
  )

set(TEMP_DIR ${ImageMath_BINARY_DIR}/Testing/Data)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Data )
set(MultiOtsuOutput ${TEMP_DIR}/MultiOtsuOutput.nrrd)


add_executable(ImageMathTest ImageMathTest.cxx)
target_link_libraries(ImageMathTest ImageMathLib)

#Normal use of MultiOtsu functionality
ExternalData_Add_Test(ImageMathData NAME MultiOtsuNormalRun COMMAND $<TARGET_FILE:ImageMathTest>
#add_test(NAME MultiOtsuNormalRun COMMAND $<TARGET_FILE:ImageMathTest>
  --compare
  DATA{Data/MultiOtsuBaseline.nrrd}
  #${CMAKE_CURRENT_SOURCE_DIR}/Data/MultiOtsuBaseline.nrrd}    
  ${MultiOtsuOutput}
  ModuleEntryPoint
  DATA{Data/MultiOtsuInput.nrrd}
  -outfile ${MultiOtsuOutput}
  -otsuMultipleThresholds
  -otsuPara
  3,0,128
)

# Add a build target to populate the real data.
ExternalData_Add_Target(ImageMathData)
