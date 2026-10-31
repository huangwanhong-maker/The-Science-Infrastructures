# The architecture domain's C++: buildings (ADR 0060 §2.8, ADR 0062).
# Header-only; each module a target architecture::<module>:
#   space           a building as a complex of rooms: its rooms the 3-cells, walls, floors and facades the 2-cells,
#                   what each room carries; doors and stairs opened, rooms merged and parted, the rooms reached;
#   environment     the heat of its rooms as the Hodge Laplacian Δ₃ with the fabric in its stars; the Sun's sol-air
#                   gains on its facades and roof; its thermal modes; each hour propagated exactly.
#
# The domain stands on the foundations alone (module.json): calx (linalg), morphe (topology, formal) and physis
# (discretisation). Inside The Science Infrastructures they are the sibling foundations, added here when no larger
# project has added them.

cmake_minimum_required(VERSION 3.20)
project(architecture LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(ARCHITECTURE_TOP_LEVEL ON)
else()
  set(ARCHITECTURE_TOP_LEVEL OFF)
endif()
option(ARCHITECTURE_CHECKS "Build the architecture domain's checks" ${ARCHITECTURE_TOP_LEVEL})

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
    message(FATAL_ERROR "the architecture domain needs the foundations: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_ROOT}/foundations/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::discretisation)
  add_subdirectory(${SCIENCE_ROOT}/foundations/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(ARCHITECTURE_CHECKS)
  enable_testing()
endif()

function(architecture_module name)
  cmake_parse_arguments(ARG "" "" "LIBRARIES" ${ARGN})
  add_library(architecture_${name} INTERFACE)
  add_library(architecture::${name} ALIAS architecture_${name})
  target_include_directories(architecture_${name} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${name}/include)
  target_compile_features(architecture_${name} INTERFACE cxx_std_20)
  target_link_libraries(architecture_${name} INTERFACE ${ARG_LIBRARIES})
  if(ARCHITECTURE_CHECKS)
    calx_add_check(check_architecture_${name} SOURCES ${name}/tests/check_${name}.cpp LIBRARIES architecture::${name})
  endif()
endfunction()

architecture_module(space LIBRARIES morphe::topology morphe::formal)
architecture_module(environment LIBRARIES architecture::space physis::discretisation calx::linalg)
