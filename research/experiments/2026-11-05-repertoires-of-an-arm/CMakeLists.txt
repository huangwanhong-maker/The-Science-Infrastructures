# Repertoires of an arm: quality-diversity on heuresis's planar arm of twelve joints.
#
# A project of its own, built against the foundations by relative path (calx, morphe, physis, tyche, heuresis),
# each added once with an explicit binary directory, as foundations/heuresis/CMakeLists.txt adds its siblings.
#
#   cmake -S . -B ~/.cache/science-build/exp-arm -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-arm -j4
#   ~/.cache/science-build/exp-arm/repertoires results

cmake_minimum_required(VERSION 3.20)
project(experiment_repertoires_of_an_arm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()
find_package(Threads REQUIRED)

set(FOUNDATIONS ${CMAKE_CURRENT_SOURCE_DIR}/../../../foundations)
if(NOT TARGET calx::linalg)
  add_subdirectory(${FOUNDATIONS}/calx ${CMAKE_BINARY_DIR}/calx)
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${FOUNDATIONS}/morphe ${CMAKE_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::units)
  add_subdirectory(${FOUNDATIONS}/physis ${CMAKE_BINARY_DIR}/physis)
endif()
if(NOT TARGET tyche::probability)
  add_subdirectory(${FOUNDATIONS}/tyche ${CMAKE_BINARY_DIR}/tyche)
endif()
if(NOT TARGET heuresis::search)
  add_subdirectory(${FOUNDATIONS}/heuresis ${CMAKE_BINARY_DIR}/heuresis)
endif()

add_executable(repertoires repertoires.cpp)
target_link_libraries(repertoires PRIVATE heuresis::diversity heuresis::strategies heuresis::benchmarks tyche::descriptive Threads::Threads)
target_compile_options(repertoires PRIVATE -Wall -Wextra -Wpedantic)
