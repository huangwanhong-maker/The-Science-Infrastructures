# A district's economy without a model of its flow: study (1) of the shared problem
# (../2026-11-05-a-district-economy-the-problem), on the search foundation heuresis.
#
# A standalone project: the foundations (calx, morphe, physis, tyche, heuresis) and the domains (architecture, then
# urban science) are added from the infrastructure by relative path, each once with an explicit binary directory, as
# the designing study's CMakeLists.txt does; the shared problem's header and the API's JSON are header-only.
#
#   cmake -S . -B ~/.cache/science-build/exp-district-economy-without -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-district-economy-without -j4
#   (from the repository's root; P = research/experiments/2026-11-05-a-district-economy-the-problem/results)
#   ~/.cache/science-build/exp-district-economy-without/without_a_model $P/problem-inputs.json $P/problem.json $P/sensors.json research/experiments/2026-11-05-a-district-economy-without-a-model/results

cmake_minimum_required(VERSION 3.20)
project(experiment_district_economy_without_a_model LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()
find_package(Threads REQUIRED)

set(SCIENCE ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
set(PROBLEM ${CMAKE_CURRENT_SOURCE_DIR}/../2026-11-05-a-district-economy-the-problem)
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
if(NOT TARGET heuresis::search)
  add_subdirectory(${SCIENCE}/foundations/heuresis ${CMAKE_BINARY_DIR}/heuresis)
endif()
if(NOT TARGET architecture::space)
  add_subdirectory(${SCIENCE}/domains/architecture ${CMAKE_BINARY_DIR}/architecture)
endif()
if(NOT TARGET urban::landscape)
  add_subdirectory(${SCIENCE}/domains/urban-science ${CMAKE_BINARY_DIR}/urban-science)
endif()

add_executable(without_a_model without_a_model.cpp)
target_include_directories(without_a_model PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${PROBLEM} ${SCIENCE}/api/cpp)
target_link_libraries(without_a_model PRIVATE urban::events urban::landscape urban::affinity architecture::space
                                              physis::statistical_mechanics physis::computational_geometry
                                              tyche::spatial tyche::regression
                                              heuresis::genetic heuresis::strategies heuresis::models heuresis::diversity heuresis::multiobjective
                                              heuresis::variation heuresis::search calx::linalg calx::random calx::scheduler Threads::Threads)
target_compile_options(without_a_model PRIVATE -Wall -Wextra)
