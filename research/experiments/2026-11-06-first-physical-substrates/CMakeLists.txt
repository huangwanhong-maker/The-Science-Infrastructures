# The first physical substrates: examples of the fluids, chemical and molecular families' simulators (ADR 0076),
# run and written to results/ for the figures of the gallery's section on physical computing.
#
# A standalone project: the foundations and the physical-computing domain are added from the infrastructure by
# relative path, each once with an explicit binary directory.
#
#   cmake -S . -B ~/.cache/science-build/exp-physical-substrates -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-physical-substrates -j4
#   ~/.cache/science-build/exp-physical-substrates/substrates <results directory>

cmake_minimum_required(VERSION 3.20)
project(experiment_first_physical_substrates LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(SCIENCE ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
if(NOT TARGET calx::verification)
  add_subdirectory(${SCIENCE}/foundations/calx ${CMAKE_BINARY_DIR}/calx)
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE}/foundations/morphe ${CMAKE_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::units)
  add_subdirectory(${SCIENCE}/foundations/physis ${CMAKE_BINARY_DIR}/physis)
endif()
if(NOT TARGET tyche::regression)
  add_subdirectory(${SCIENCE}/foundations/tyche ${CMAKE_BINARY_DIR}/tyche)
endif()
add_subdirectory(${SCIENCE}/domains/physical-computing ${CMAKE_BINARY_DIR}/physical-computing)

add_executable(substrates substrates.cpp)
target_link_libraries(substrates PRIVATE physical_computing::fluids physical_computing::chemical physical_computing::molecular)
target_compile_options(substrates PRIVATE -Wall -Wextra)
