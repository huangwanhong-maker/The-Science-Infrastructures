# Quantum search, as the search foundation heuresis states it: Grover's law and its exact emulation, Grover
# adaptive search in oracle calls, and the quantum-inspired EA as a product-Bernoulli search distribution.
#
# A standalone project: the foundations are added from the infrastructure by relative path, each once (as
# foundations/heuresis/CMakeLists.txt adds its siblings).

cmake_minimum_required(VERSION 3.20)
project(quantum_search LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

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
if(NOT TARGET heuresis::search)
  add_subdirectory(${SCIENCE}/foundations/heuresis ${CMAKE_BINARY_DIR}/heuresis)
endif()

add_executable(quantum quantum.cpp)
target_compile_options(quantum PRIVATE -Wall -Wextra)
target_link_libraries(quantum PRIVATE heuresis::quantum heuresis::benchmarks)
