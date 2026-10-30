# The oceanography domain's C++ (ADR 0059): its first module, water, the ocean's water as a body in it meets
# it, its density and near-surface current presumed and realised as fields. Header-only; each module a target
# oceanography::<module>. The domain's provisional JavaScript (js/) is checked by tests/run-tests.mjs with
# these checks.
#
# The domain may use geodesy (module.json); water uses only the foundations: physis (stationary random fields,
# physis::analysis; fields, physis::fields) and calx (the report of its checks). Inside The Science
# Infrastructures they are the sibling foundations, added here when no larger project has added them.

cmake_minimum_required(VERSION 3.20)
project(oceanography LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(OCEANOGRAPHY_TOP_LEVEL ON)
else()
  set(OCEANOGRAPHY_TOP_LEVEL OFF)
endif()
option(OCEANOGRAPHY_CHECKS "Build the oceanography domain's checks" ${OCEANOGRAPHY_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(SCIENCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)
if(NOT TARGET calx::verification)
  if(EXISTS ${SCIENCE_ROOT}/foundations/calx/CMakeLists.txt)
    add_subdirectory(${SCIENCE_ROOT}/foundations/calx ${CMAKE_CURRENT_BINARY_DIR}/calx)
  else()
    message(FATAL_ERROR "the oceanography domain needs calx and physis: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_ROOT}/foundations/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::analysis)
  add_subdirectory(${SCIENCE_ROOT}/foundations/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(OCEANOGRAPHY_CHECKS)
  enable_testing()
endif()

# water: the ocean's water as a body meets it (presumed.hpp).
add_library(oceanography_water INTERFACE)
add_library(oceanography::water ALIAS oceanography_water)
target_include_directories(oceanography_water INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/water/include)
target_compile_features(oceanography_water INTERFACE cxx_std_20)
target_link_libraries(oceanography_water INTERFACE physis::analysis physis::fields)

if(OCEANOGRAPHY_CHECKS)
  calx_add_check(check_oceanography_water SOURCES water/tests/check_water.cpp LIBRARIES oceanography::water)
endif()
