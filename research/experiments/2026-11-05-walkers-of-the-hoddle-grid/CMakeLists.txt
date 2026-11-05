# The walkers of the Hoddle Grid: urban science's walkers as agents (domains/urban-science/movement) on the district of
# the shared problem "A district's economy", against the flow study's logit flow, with a street closed, a festival's
# crowd and the street's regime decided online.
#
# A standalone project: the foundations (calx, morphe, physis, tyche, nomos, heuresis) and the domains (architecture,
# urban science, economics) are added from the infrastructure by relative path, each once with an explicit binary
# directory; the shared problem's directory is added for its INTERFACE target `district_economy`, and the flow study's
# headers (flow.hpp) are included from its directory: its calibrated trips and its logit flow are this study's
# macroscopic reference.
#
#   cmake -S . -B ~/.cache/science-build/exp-walkers -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-walkers -j2
#   (from the repository's root)
#   ~/.cache/science-build/exp-walkers/walkers <problem-inputs.json> <problem.json> <sensors.json> <results directory> [--threads N] [--task week|closure|festival|regime|all]

cmake_minimum_required(VERSION 3.20)
project(experiment_walkers_of_the_hoddle_grid LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()
find_package(Threads REQUIRED)

set(SCIENCE ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
set(PROBLEM ${CMAKE_CURRENT_SOURCE_DIR}/../2026-11-05-a-district-economy-the-problem)
set(FLOW ${CMAKE_CURRENT_SOURCE_DIR}/../2026-11-05-a-district-economy-as-a-flow)
if(NOT TARGET calx::linalg)
  add_subdirectory(${SCIENCE}/foundations/calx ${CMAKE_BINARY_DIR}/calx)
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE}/foundations/morphe ${CMAKE_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::units)
  add_subdirectory(${SCIENCE}/foundations/physis ${CMAKE_BINARY_DIR}/physis)
endif()
if(NOT TARGET tyche::probability)
  add_subdirectory(${SCIENCE}/foundations/tyche ${CMAKE_BINARY_DIR}/tyche)
endif()
if(NOT TARGET nomos::measurement)
  add_subdirectory(${SCIENCE}/foundations/nomos ${CMAKE_BINARY_DIR}/nomos)
endif()
if(NOT TARGET heuresis::search)
  add_subdirectory(${SCIENCE}/foundations/heuresis ${CMAKE_BINARY_DIR}/heuresis)
endif()
if(NOT TARGET architecture::space)
  add_subdirectory(${SCIENCE}/domains/architecture ${CMAKE_BINARY_DIR}/architecture)
endif()
if(NOT TARGET urban::landscape)
  add_subdirectory(${SCIENCE}/domains/urban-science ${CMAKE_BINARY_DIR}/urban-science)
endif()
if(NOT TARGET economics::circulation)
  add_subdirectory(${SCIENCE}/domains/economics ${CMAKE_BINARY_DIR}/economics)
endif()
if(NOT TARGET district_economy)
  add_subdirectory(${PROBLEM} ${CMAKE_BINARY_DIR}/problem EXCLUDE_FROM_ALL)
endif()

add_executable(walkers walkers.cpp)
target_include_directories(walkers PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${FLOW})
target_link_libraries(walkers PRIVATE district_economy urban::movement urban::flow urban::rhythm economics::circulation physis::discretisation physis::computational_geometry
                                      tyche::regression tyche::decision heuresis::genetic heuresis::strategies heuresis::search heuresis::variation calx::linalg calx::random
                                      calx::scheduler Threads::Threads)
target_compile_options(walkers PRIVATE -Wall -Wextra)
