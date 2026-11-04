# Designing a district: the redesign of a few blocks of Melbourne's Hoddle Grid as a heterogeneous, many-objective
# search, on the search foundation heuresis and the architecture and urban-science domains.
#
# A standalone project: the foundations (calx, morphe, physis, tyche, heuresis) and the domains (architecture, then
# urban science, which stands on it) are added from the infrastructure by relative path, each once with an explicit
# binary directory, as foundations/heuresis/CMakeLists.txt adds its siblings; the API's JSON is header-only.
#
#   cmake -S . -B ~/.cache/science-build/exp-designing-a-district -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-designing-a-district -j4
#   ~/.cache/science-build/exp-designing-a-district/designing_a_district results/district-inputs.json results

cmake_minimum_required(VERSION 3.20)
project(experiment_designing_a_district LANGUAGES CXX)

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
if(NOT TARGET heuresis::search)
  add_subdirectory(${SCIENCE}/foundations/heuresis ${CMAKE_BINARY_DIR}/heuresis)
endif()
if(NOT TARGET architecture::space)
  add_subdirectory(${SCIENCE}/domains/architecture ${CMAKE_BINARY_DIR}/architecture)
endif()
if(NOT TARGET urban::landscape)
  add_subdirectory(${SCIENCE}/domains/urban-science ${CMAKE_BINARY_DIR}/urban-science)
endif()

add_executable(designing_a_district designing_a_district.cpp)
target_include_directories(designing_a_district PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${SCIENCE}/api/cpp)
target_link_libraries(designing_a_district PRIVATE heuresis::multiobjective heuresis::diversity heuresis::variation heuresis::search
                                                   urban::landscape architecture::space architecture::environment
                                                   physis::statistical_mechanics physis::computational_geometry calx::linalg calx::scheduler Threads::Threads)
target_compile_options(designing_a_district PRIVATE -Wall -Wextra)
