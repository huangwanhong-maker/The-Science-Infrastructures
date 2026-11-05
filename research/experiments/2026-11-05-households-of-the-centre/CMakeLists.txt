# The households of the centre (ADR 0075): a case study of urban-science's housing on the City of Melbourne.
#
# A standalone project: the foundations (calx, morphe, physis, tyche, nomos) and the domains (architecture, urban
# science) are added from the infrastructure by relative path, each once with an explicit binary directory; the API's
# JSON (api/cpp) reads the data library's product.
#
#   cmake -S . -B ~/.cache/science-build/exp-households -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-households -j2
#   (from the repository's root, after python data/products/housing_records.py)
#   ~/.cache/science-build/exp-households/households data/cache/products/housing-melbourne.json research/experiments/2026-11-05-households-of-the-centre/results

cmake_minimum_required(VERSION 3.20)
project(experiment_households_of_the_centre LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()
find_package(Threads REQUIRED)

set(SCIENCE ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
set(URBAN_SCIENCE_CHECKS OFF CACHE BOOL "" FORCE)
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
if(NOT TARGET urban::housing)
  add_subdirectory(${SCIENCE}/domains/urban-science ${CMAKE_BINARY_DIR}/urban-science)
endif()

add_executable(households households.cpp)
target_include_directories(households PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${SCIENCE}/api/cpp)
target_link_libraries(households PRIVATE urban::housing tyche::decision calx::signal calx::random Threads::Threads)
target_compile_options(households PRIVATE -Wall -Wextra)
