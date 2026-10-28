# The metrology domain's C++: time scales and the calendar (time), the Earth orientation
# parameters (earth_orientation), and the orientation of the Earth in space with the geocentric
# frames (frames). Header-only; each module a target metrology::<module>.
#
# The domain uses physis (rotations, rigid motions and the frame graph, physis::geometry) and calx
# (the report of its checks). Inside The Science Infrastructures they are the sibling foundations,
# added here when no larger project has added them.

cmake_minimum_required(VERSION 3.20)
project(metrology LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(METROLOGY_TOP_LEVEL ON)
else()
  set(METROLOGY_TOP_LEVEL OFF)
endif()
option(METROLOGY_CHECKS "Build the metrology domain's checks" ${METROLOGY_TOP_LEVEL})

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
    message(FATAL_ERROR "the metrology domain needs calx and physis: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::geometry)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(METROLOGY_CHECKS)
  enable_testing()
endif()

# earth_orientation: the IERS's Earth orientation parameters, read and interpolated.
add_library(metrology_earth_orientation INTERFACE)
add_library(metrology::earth_orientation ALIAS metrology_earth_orientation)
target_include_directories(metrology_earth_orientation INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/earth_orientation/include)
target_compile_features(metrology_earth_orientation INTERFACE cxx_std_20)

# time: Julian dates, the calendar, leap seconds, the time scales and their conversions.
add_library(metrology_time INTERFACE)
add_library(metrology::time ALIAS metrology_time)
target_include_directories(metrology_time INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/time/include)
target_link_libraries(metrology_time INTERFACE metrology::earth_orientation)

# frames: IAU 2006/2000A precession-nutation, Earth rotation, polar motion; the geocentric frames.
add_library(metrology_frames INTERFACE)
add_library(metrology::frames ALIAS metrology_frames)
target_include_directories(metrology_frames INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/frames/include)
target_link_libraries(metrology_frames INTERFACE metrology::time physis::geometry)

if(METROLOGY_CHECKS)
  set(METROLOGY_FIXTURES ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures)
  calx_add_check(check_metrology_time SOURCES time/tests/check_time.cpp LIBRARIES metrology::time)
  target_compile_definitions(check_metrology_time PRIVATE METROLOGY_FIXTURES="${METROLOGY_FIXTURES}")
  calx_add_check(check_metrology_frames SOURCES frames/tests/check_frames.cpp LIBRARIES metrology::frames)
  target_compile_definitions(check_metrology_frames PRIVATE METROLOGY_FIXTURES="${METROLOGY_FIXTURES}")
endif()
