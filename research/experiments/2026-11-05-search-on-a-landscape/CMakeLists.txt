# Search on a landscape: how six algorithms of heuresis move over three two-dimensional functions.
#
# A standalone project on the foundations under ../../../foundations, as heuresis adds its siblings: each library
# is added once, when no larger project has added it. Build (in WSL, from this directory):
#
#   cmake -S . -B ~/.cache/science-build/exp-search-on-a-landscape -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-search-on-a-landscape -j4
#   ~/.cache/science-build/exp-search-on-a-landscape/landscape results

cmake_minimum_required(VERSION 3.20)
project(search_on_a_landscape LANGUAGES CXX)

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

# The commit of heuresis, written into the results.
execute_process(COMMAND git -C "${FOUNDATIONS}/heuresis" rev-parse --short HEAD OUTPUT_VARIABLE HEURESIS_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

add_executable(landscape landscape.cpp)
target_compile_options(landscape PRIVATE -Wall -Wextra -Wpedantic)
if(HEURESIS_COMMIT)
  target_compile_definitions(landscape PRIVATE HEURESIS_COMMIT="${HEURESIS_COMMIT}")
endif()
target_link_libraries(landscape PRIVATE
  heuresis::search heuresis::benchmarks heuresis::strategies heuresis::swarms heuresis::distributions
  heuresis::models heuresis::analysis Threads::Threads)
