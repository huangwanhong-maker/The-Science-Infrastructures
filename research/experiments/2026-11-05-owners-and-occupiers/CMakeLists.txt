# Owners and occupiers of the Hoddle Grid's buildings: a case study of the urban-science domain's `occupancy` (ADR 0071,
# ADR 0075). A standalone project: the foundations (calx, morphe, physis, tyche, nomos) and the domains (architecture,
# urban science) are added from the infrastructure by relative path, each once with an explicit binary directory; the
# API's C++ (science/json.hpp) is read for the record's JSON.
#
#   cmake -S . -B ~/.cache/science-build/exp-owners-occupiers -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-owners-occupiers -j4
#   (from the repository's root)
#   ~/.cache/science-build/exp-owners-occupiers/owners_and_occupiers data/cache/products/occupancy-hoddle-grid.json \
#       data/cache/products/city-events-melbourne.json research/experiments/2026-11-05-owners-and-occupiers/results

cmake_minimum_required(VERSION 3.20)
project(experiment_owners_and_occupiers LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()
find_package(Threads REQUIRED)

set(SCIENCE ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
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
if(NOT TARGET nomos::agents)
  add_subdirectory(${SCIENCE}/foundations/nomos ${CMAKE_BINARY_DIR}/nomos)
endif()
if(NOT TARGET architecture::space)
  add_subdirectory(${SCIENCE}/domains/architecture ${CMAKE_BINARY_DIR}/architecture)
endif()
if(NOT TARGET urban::occupancy)
  add_subdirectory(${SCIENCE}/domains/urban-science ${CMAKE_BINARY_DIR}/urban-science)
endif()

add_executable(owners_and_occupiers owners_and_occupiers.cpp)
target_include_directories(owners_and_occupiers PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${SCIENCE}/api/cpp)
target_link_libraries(owners_and_occupiers PRIVATE urban::occupancy calx::scheduler Threads::Threads)
target_compile_options(owners_and_occupiers PRIVATE -Wall -Wextra)
