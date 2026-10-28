# The geophysics domain's C++: the geomagnetic field (geomagnetism: the World Magnetic Model 2025,
# its elements and secular variation at a geodetic place, the field on the Earth-fixed axes, the
# centred dipole and its L, the model's validity and error model). Header-only; each module a target
# geophysics::<module>. The domain's plates, mechanisms, crust and aftershocks are still JavaScript
# (js/). The checks read WMM2025's coefficients from the data library's cache (noaa-wmm/WMM.COF, or
# under SCIENCE_CACHE when set) and skip them, saying so, when they are missing.
#
# The domain uses geodesy (the figure of the Earth: WGS 84, geodetic coordinates, local frames),
# metrology (the calendar, for decimal years), physis (vectors and rotations, physis::geometry) and
# calx (special functions, random streams and the report of its checks). Inside The Science
# Infrastructures they are the sibling foundations and domains, added here when no larger project
# has added them.

cmake_minimum_required(VERSION 3.20)
project(geophysics LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(GEOPHYSICS_TOP_LEVEL ON)
else()
  set(GEOPHYSICS_TOP_LEVEL OFF)
endif()
option(GEOPHYSICS_CHECKS "Build the geophysics domain's checks" ${GEOPHYSICS_TOP_LEVEL})

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
    message(FATAL_ERROR "the geophysics domain needs calx and physis: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_ROOT}/foundations/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::geometry)
  add_subdirectory(${SCIENCE_ROOT}/foundations/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
# The domains beneath, in their order, when no larger project has added them.
if(NOT TARGET metrology::time)
  add_subdirectory(${SCIENCE_ROOT}/domains/metrology ${CMAKE_CURRENT_BINARY_DIR}/metrology)
endif()
if(NOT TARGET geodesy::figure)
  add_subdirectory(${SCIENCE_ROOT}/domains/geodesy ${CMAKE_CURRENT_BINARY_DIR}/geodesy)
endif()
if(GEOPHYSICS_CHECKS)
  enable_testing()
endif()

# geomagnetism: the World Magnetic Model (WMM2025): coefficients, the field and its secular variation,
# the centred dipole and L, validity and the error model (the Schmidt functions are calx's).
add_library(geophysics_geomagnetism INTERFACE)
add_library(geophysics::geomagnetism ALIAS geophysics_geomagnetism)
target_include_directories(geophysics_geomagnetism INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/geomagnetism/include)
target_compile_features(geophysics_geomagnetism INTERFACE cxx_std_20)
target_link_libraries(geophysics_geomagnetism INTERFACE geodesy::figure metrology::time physis::geometry calx::special_functions)

if(GEOPHYSICS_CHECKS)
  calx_add_check(check_geophysics_wmm SOURCES geomagnetism/tests/check_wmm.cpp LIBRARIES geophysics::geomagnetism calx::random)
  target_compile_definitions(check_geophysics_wmm PRIVATE GEOPHYSICS_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures"
                             GEOPHYSICS_ROOT="${CMAKE_CURRENT_SOURCE_DIR}")
endif()
