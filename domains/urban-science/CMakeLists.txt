# The urban-science domain's C++: streets as open systems, after the user's Street Oceanography (ADR 0060).
# Header-only; each module a target urban::<module>:
#   street          the street system: holdings, uses and declared hours; terrain, its contractions and its
#                   dominance order; the length of the terrain and the compression ratio; the departure from
#                   mixing and its scales; the effective terrain in time; the visits of kinds of people;
#   rhythm          calendars as angles and the lattice of their combinations; records of counts in UTC and
#                   civil time; spectra against the integration null, lines read on the lattice, harmonic
#                   constants, selectivity, the civil share of the clock, bicoherence, the paper's Figure 1;
#   response        susceptibility to irregular and dated forcings; the comparison across a repeating alteration;
#   affinity        co-location of uses against free and terrain-stratified nulls; catalysis in space and time
#                   (the bridal store and the ring store); the spacing of clusters;
#   flow            the crowd as a kinematic wave under a rule that opens a carriageway; person-time; movement
#                   among places, its transit times and its circulation;
#   classification  the plane of compression and selectivity, and Clark and Evans's test.
#
# The domain stands on the foundations alone (module.json): calx (signal, linalg, random), morphe (topology),
# physis (continuum, statistical_mechanics, chemistry, discretisation) and tyche (time_series, regression,
# spatial). Inside The Science Infrastructures they are the sibling foundations, added here when no larger
# project has added them.

cmake_minimum_required(VERSION 3.20)
project(urban_science LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(URBAN_SCIENCE_TOP_LEVEL ON)
else()
  set(URBAN_SCIENCE_TOP_LEVEL OFF)
endif()
option(URBAN_SCIENCE_CHECKS "Build the urban-science domain's checks" ${URBAN_SCIENCE_TOP_LEVEL})

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
    message(FATAL_ERROR "the urban-science domain needs the foundations: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_ROOT}/foundations/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::statistical_mechanics)
  add_subdirectory(${SCIENCE_ROOT}/foundations/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(NOT TARGET tyche::time_series)
  add_subdirectory(${SCIENCE_ROOT}/foundations/tyche ${CMAKE_CURRENT_BINARY_DIR}/tyche)
endif()
if(URBAN_SCIENCE_CHECKS)
  enable_testing()
endif()

function(urban_module name)
  cmake_parse_arguments(ARG "" "" "LIBRARIES" ${ARGN})
  add_library(urban_${name} INTERFACE)
  add_library(urban::${name} ALIAS urban_${name})
  target_include_directories(urban_${name} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${name}/include)
  target_compile_features(urban_${name} INTERFACE cxx_std_20)
  target_link_libraries(urban_${name} INTERFACE ${ARG_LIBRARIES})
  if(URBAN_SCIENCE_CHECKS)
    calx_add_check(check_urban_${name} SOURCES ${name}/tests/check_${name}.cpp LIBRARIES urban::${name})
  endif()
endfunction()

urban_module(street LIBRARIES physis::statistical_mechanics calx::linalg calx::random)
urban_module(rhythm LIBRARIES tyche::time_series calx::signal calx::random)
urban_module(response LIBRARIES urban::rhythm tyche::regression calx::random)
urban_module(affinity LIBRARIES tyche::spatial tyche::time_series physis::statistical_mechanics physis::chemistry calx::linalg)
urban_module(flow LIBRARIES urban::rhythm tyche::probability physis::continuum physis::statistical_mechanics physis::discretisation morphe::topology calx::linalg)
urban_module(classification LIBRARIES tyche::spatial)
