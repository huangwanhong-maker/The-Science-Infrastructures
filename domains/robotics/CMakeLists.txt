# The robotics domain's C++: systems that sense and act, their bodies, controllers, plans and sensors (ADR
# 0054). Header-only; each module a target robotics::<module>:
#   vehicles  a vehicle on a sphere at two fidelities, kinematic and dynamic, and the translation between them;
#   control   following waypoints;
#   planning  surveys: a lawnmower pattern, rings of candidates, the greedy choice of the next place;
#   sensing   mounts and white noise;
#   legged    a rabbit: a hopping robot's body, motors, tendon and feet in its plane, its controller (standing,
#             hopping in place and forward, stopping), and its simulation by three methods of one motion (ADR 0058).
#
# The domain stands on the foundations alone (module.json): physis (rotations, physis::geometry; articulated
# bodies and contact, physis::mechanics) and calx
# (RATTLE, calx::integrators; counter-based streams, calx::random; the report of its checks). Inside The
# Science Infrastructures they are the sibling foundations, added here when no larger project has added them.

cmake_minimum_required(VERSION 3.20)
project(robotics LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(ROBOTICS_TOP_LEVEL ON)
else()
  set(ROBOTICS_TOP_LEVEL OFF)
endif()
option(ROBOTICS_CHECKS "Build the robotics domain's checks" ${ROBOTICS_TOP_LEVEL})

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
    message(FATAL_ERROR "the robotics domain needs calx and physis: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_ROOT}/foundations/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::geometry)
  add_subdirectory(${SCIENCE_ROOT}/foundations/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(NOT TARGET tyche::spatial)
  add_subdirectory(${SCIENCE_ROOT}/foundations/tyche ${CMAKE_CURRENT_BINARY_DIR}/tyche)
endif()
if(ROBOTICS_CHECKS)
  enable_testing()
endif()

function(robotics_module name)
  cmake_parse_arguments(ARG "" "" "LIBRARIES" ${ARGN})
  add_library(robotics_${name} INTERFACE)
  add_library(robotics::${name} ALIAS robotics_${name})
  target_include_directories(robotics_${name} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${name}/include)
  target_compile_features(robotics_${name} INTERFACE cxx_std_20)
  target_link_libraries(robotics_${name} INTERFACE ${ARG_LIBRARIES})
  if(ROBOTICS_CHECKS)
    calx_add_check(check_robotics_${name} SOURCES ${name}/tests/check_${name}.cpp LIBRARIES robotics::${name})
  endif()
endfunction()

robotics_module(vehicles LIBRARIES physis::geometry calx::integrators)
robotics_module(actuation)
robotics_module(control LIBRARIES robotics::vehicles)
robotics_module(planning LIBRARIES robotics::vehicles)
robotics_module(sensing LIBRARIES robotics::vehicles calx::random)
robotics_module(legged LIBRARIES robotics::actuation physis::mechanics calx::integrators calx::differentiation calx::linalg calx::random tyche::spatial physis::analysis)
robotics_module(swimming LIBRARIES robotics::actuation physis::mechanics calx::integrators calx::differentiation calx::linalg)
robotics_module(quadruped LIBRARIES robotics::actuation robotics::legged physis::mechanics calx::integrators calx::differentiation calx::linalg calx::random)
if(ROBOTICS_CHECKS)
  target_link_libraries(check_robotics_control PRIVATE robotics::planning)
endif()
