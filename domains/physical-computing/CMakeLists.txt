# The physical-computing domain's C++: computation by the dynamics of matter (ADR 0076). Header-only; each
# module a target physical_computing::<module>, in the pipeline's order:
#   tasks      what is computed, substrate-neutral (tasks/);
#   schedule   when inputs are applied and states read (schedule/);
#   computer   the commuting square, its tolerance, the cost of a verified result and the materials' marks (computer/);
#   runtime    the backend interface, provenance, and programs driven on backends (runtime/);
#   fluids     the fluids family's simulators (substrates/fluids/);
#   chemical   the chemical family's: the BZ medium, its arrays, reaction networks as programs, electrodes (substrates/chemical/);
#   molecular  the molecular family's: strand displacement, seesaw gates, the PEN predator–prey (substrates/molecular/).
# The intermediate representation (ir/) is optional and future, and has no code.
#
# The domain stands on the foundations alone (module.json): calx, morphe, physis and tyche, the sibling
# foundations inside The Science Infrastructures, added here when no larger project has added them.

cmake_minimum_required(VERSION 3.20)
project(physical_computing LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(PHYSICAL_COMPUTING_TOP_LEVEL ON)
else()
  set(PHYSICAL_COMPUTING_TOP_LEVEL OFF)
endif()
option(PHYSICAL_COMPUTING_CHECKS "Build the physical-computing domain's checks" ${PHYSICAL_COMPUTING_TOP_LEVEL})

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
    message(FATAL_ERROR "the physical-computing domain needs its foundations: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_ROOT}/foundations/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::dynamics)
  add_subdirectory(${SCIENCE_ROOT}/foundations/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(NOT TARGET tyche::regression)
  add_subdirectory(${SCIENCE_ROOT}/foundations/tyche ${CMAKE_CURRENT_BINARY_DIR}/tyche)
endif()
if(PHYSICAL_COMPUTING_CHECKS)
  enable_testing()
endif()

# A module: its headers under <dir>/include, its checks <dir>/tests/check_<name>.cpp for each name in CHECKS.
function(physical_computing_module name)
  cmake_parse_arguments(ARG "" "DIR" "LIBRARIES;CHECKS" ${ARGN})
  if(NOT ARG_DIR)
    set(ARG_DIR ${name})
  endif()
  add_library(physical_computing_${name} INTERFACE)
  add_library(physical_computing::${name} ALIAS physical_computing_${name})
  target_include_directories(physical_computing_${name} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_DIR}/include)
  target_compile_features(physical_computing_${name} INTERFACE cxx_std_20)
  target_link_libraries(physical_computing_${name} INTERFACE ${ARG_LIBRARIES})
  if(PHYSICAL_COMPUTING_CHECKS)
    foreach(check ${ARG_CHECKS})
      calx_add_check(check_physical_computing_${check} SOURCES ${ARG_DIR}/tests/check_${check}.cpp LIBRARIES physical_computing::${name})
    endforeach()
  endif()
endfunction()

physical_computing_module(tasks)
physical_computing_module(schedule)
physical_computing_module(computer LIBRARIES physical_computing::tasks physis::statistical_mechanics CHECKS computer)
physical_computing_module(runtime LIBRARIES physical_computing::computer physical_computing::schedule CHECKS runtime)
physical_computing_module(fluids DIR substrates/fluids
  LIBRARIES physical_computing::runtime physis::continuum physis::dynamics physis::discretisation physis::statistical_mechanics calx::random calx::optimisation calx::solvers
  CHECKS tank convection_cell rotating_tank faraday walkers capillary foam bubble_logic phillips)
physical_computing_module(chemical DIR substrates/chemical
  LIBRARIES physical_computing::runtime physis::chemistry physis::dynamics calx::solvers
  CHECKS oregonator bz_array programs electrochemical)
physical_computing_module(molecular DIR substrates/molecular
  LIBRARIES physical_computing::runtime physis::chemistry physis::dynamics
  CHECKS strand_displacement seesaw pen)
