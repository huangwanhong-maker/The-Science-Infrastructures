# Benchmarking continuous search: ten algorithms of heuresis on eight shifted and rotated functions in 2, 5, 10 and
# 20 dimensions, in the manner of COCO.
#
# A standalone project on the foundations under ../../../foundations, as heuresis adds its siblings: each library
# is added once, when no larger project has added it. Build (in WSL, from this directory):
#
#   cmake -S . -B ~/.cache/science-build/exp-benchmarking-continuous-search -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-benchmarking-continuous-search -j4
#   ~/.cache/science-build/exp-benchmarking-continuous-search/benchmark results

cmake_minimum_required(VERSION 3.20)
project(benchmarking_continuous_search LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
find_package(Threads REQUIRED)

set(FOUNDATIONS "${CMAKE_CURRENT_SOURCE_DIR}/../../../foundations")
if(NOT TARGET calx::linalg)
  add_subdirectory("${FOUNDATIONS}/calx" "${CMAKE_BINARY_DIR}/calx")
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory("${FOUNDATIONS}/morphe" "${CMAKE_BINARY_DIR}/morphe")
endif()
if(NOT TARGET physis::units)
  add_subdirectory("${FOUNDATIONS}/physis" "${CMAKE_BINARY_DIR}/physis")
endif()
if(NOT TARGET tyche::probability)
  add_subdirectory("${FOUNDATIONS}/tyche" "${CMAKE_BINARY_DIR}/tyche")
endif()
if(NOT TARGET heuresis::search)
  add_subdirectory("${FOUNDATIONS}/heuresis" "${CMAKE_BINARY_DIR}/heuresis")
endif()

add_executable(benchmark benchmark.cpp)
target_compile_options(benchmark PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(benchmark PRIVATE
  heuresis::search heuresis::benchmarks heuresis::strategies heuresis::swarms heuresis::distributions heuresis::analysis
  Threads::Threads)
