# The astronomy domain's C++: the JPL ephemerides and binary orientation kernels (ephemerides),
# and the rotation of the bodies, the IAU models and the Moon's DE440 frames (rotation).
# Header-only; each module a target astronomy::<module>.
#
# The domain uses physis (rotations and rigid motions, physis::geometry) and calx (the report of
# its checks); it may use metrology and geodesy beneath it. Its checks read the kernels from the
# infrastructure's data cache (jpl-de440s, naif-lunar-frames) and skip what is missing.

cmake_minimum_required(VERSION 3.20)
project(astronomy LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(ASTRONOMY_TOP_LEVEL ON)
else()
  set(ASTRONOMY_TOP_LEVEL OFF)
endif()
option(ASTRONOMY_CHECKS "Build the astronomy domain's checks" ${ASTRONOMY_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(SCIENCE_FOUNDATIONS ${CMAKE_CURRENT_SOURCE_DIR}/../../foundations)
if(NOT TARGET calx::verification)
  if(EXISTS ${SCIENCE_FOUNDATIONS}/calx/CMakeLists.txt)
    add_subdirectory(${SCIENCE_FOUNDATIONS}/calx ${CMAKE_CURRENT_BINARY_DIR}/calx)
  else()
    message(FATAL_ERROR "the astronomy domain needs calx and physis: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::geometry)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(ASTRONOMY_CHECKS)
  enable_testing()
endif()

# ephemerides: DAF files, SPK ephemerides (types 2 and 3), binary PCK orientation (type 2), bodies.
add_library(astronomy_ephemerides INTERFACE)
add_library(astronomy::ephemerides ALIAS astronomy_ephemerides)
target_include_directories(astronomy_ephemerides INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/ephemerides/include)
target_compile_features(astronomy_ephemerides INTERFACE cxx_std_20)
target_link_libraries(astronomy_ephemerides INTERFACE physis::geometry)

# rotation: the IAU rotation models of the Sun, the Moon and the planets; the Moon's DE440 frames.
add_library(astronomy_rotation INTERFACE)
add_library(astronomy::rotation ALIAS astronomy_rotation)
target_include_directories(astronomy_rotation INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/rotation/include)
target_link_libraries(astronomy_rotation INTERFACE astronomy::ephemerides physis::geometry)

if(ASTRONOMY_CHECKS)
  calx_add_check(check_astronomy_ephemerides SOURCES ephemerides/tests/check_ephemerides.cpp LIBRARIES astronomy::rotation)
  target_compile_definitions(check_astronomy_ephemerides PRIVATE ASTRONOMY_ROOT="${CMAKE_CURRENT_SOURCE_DIR}")
endif()
