cmake_minimum_required(VERSION 3.5)

project(DTITensorNormalization)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

niral_add_executable(
  NAME DTITensorNormalization
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  )
