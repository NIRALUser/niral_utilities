project(DTITensorNormalization)

find_package(ITK 5.3 REQUIRED)
include(${ITK_USE_FILE})

niral_add_executable(
  NAME DTITensorNormalization
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  )
