# The planetary-science domain's C++: the planets and moons as worlds (ADR 0052). Its first module,
# magnetism: the crustal field of a body without a dynamo (the Moon's, T2015_449, with the law of its
# spectrum and the part it leaves unresolved) and the interplanetary field carried from L1 to a body
# near the Earth. Header-only; each module a target planetary_science::<module>. The checks read
# T2015_449 from the data library's cache (moon-magnetic-t2015/T2015_449.sh, or under SCIENCE_CACHE
# when set) and skip those checks, saying so, when it is missing.
#
# The domain may use astronomy, geophysics and meteorology (module.json); magnetism uses only the
# foundations: physis (potential fields and random fields, physis::analysis) and calx (the report of
# its checks). Inside The Science Infrastructures they are the sibling foundations, added here when no
# larger project has added them.

cmake_minimum_required(VERSION 3.20)
project(planetary_science LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(PLANETARY_SCIENCE_TOP_LEVEL ON)
else()
  set(PLANETARY_SCIENCE_TOP_LEVEL OFF)
endif()
option(PLANETARY_SCIENCE_CHECKS "Build the planetary-science domain's checks" ${PLANETARY_SCIENCE_TOP_LEVEL})

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
    message(FATAL_ERROR "the planetary-science domain needs calx and physis: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_ROOT}/foundations/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::analysis)
  add_subdirectory(${SCIENCE_ROOT}/foundations/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(PLANETARY_SCIENCE_CHECKS)
  enable_testing()
endif()

# magnetism: the crustal field of a body without a dynamo, and the interplanetary field carried to it.
add_library(planetary_science_magnetism INTERFACE)
add_library(planetary_science::magnetism ALIAS planetary_science_magnetism)
target_include_directories(planetary_science_magnetism INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/magnetism/include)
target_compile_features(planetary_science_magnetism INTERFACE cxx_std_20)
target_link_libraries(planetary_science_magnetism INTERFACE physis::analysis)

if(PLANETARY_SCIENCE_CHECKS)
  calx_add_check(check_planetary_science_magnetism SOURCES magnetism/tests/check_magnetism.cpp LIBRARIES planetary_science::magnetism)
  target_compile_definitions(check_planetary_science_magnetism PRIVATE PLANETARY_SCIENCE_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures"
                             PLANETARY_SCIENCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}")
endif()
