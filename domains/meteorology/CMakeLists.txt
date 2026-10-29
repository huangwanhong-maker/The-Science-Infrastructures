# The meteorology domain's C++: the U.S. Standard Atmosphere 1976 to 86 km, and its quantities as fields
# (atmosphere). Header-only;
# each module a target meteorology::<module>. The rest of the domain is provisional JavaScript (js/).
#
# The domain uses physis (fields, physis::fields) and calx (the report of its checks); it may use geodesy
# and the domains beneath it in the order of the domains. Inside The Science Infrastructures they are the
# sibling foundations, added here when no larger project has added them.

cmake_minimum_required(VERSION 3.20)
project(meteorology LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(METEOROLOGY_TOP_LEVEL ON)
else()
  set(METEOROLOGY_TOP_LEVEL OFF)
endif()
option(METEOROLOGY_CHECKS "Build the meteorology domain's checks" ${METEOROLOGY_TOP_LEVEL})

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
    message(FATAL_ERROR "the meteorology domain needs calx and physis: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::fields)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(METEOROLOGY_CHECKS)
  enable_testing()
endif()

# atmosphere: the U.S. Standard Atmosphere 1976 from -5 km' to 86 km, columns of levels through it, and
# its quantities as fields on the line of heights (field.hpp).
add_library(meteorology_atmosphere INTERFACE)
add_library(meteorology::atmosphere ALIAS meteorology_atmosphere)
target_include_directories(meteorology_atmosphere INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/atmosphere/include)
target_compile_features(meteorology_atmosphere INTERFACE cxx_std_20)
target_link_libraries(meteorology_atmosphere INTERFACE physis::fields)

if(METEOROLOGY_CHECKS)
  calx_add_check(check_meteorology_atmosphere SOURCES atmosphere/tests/check_standard.cpp LIBRARIES meteorology::atmosphere)
  target_compile_definitions(check_meteorology_atmosphere PRIVATE METEOROLOGY_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures")
endif()
