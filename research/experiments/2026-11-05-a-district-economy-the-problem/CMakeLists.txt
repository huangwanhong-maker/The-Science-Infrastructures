# A district's economy: the problem three case studies share. Its calibration (the counters' rhythm and response,
# the kinds' affinities, the district's frontages and relations) and the checks of its header, economy.hpp, which the
# studies include.
#
# A standalone project: the foundations (calx, morphe, physis, tyche, heuresis) and the domains (architecture, then
# urban science, which stands on it) are added from the infrastructure by relative path, each once with an explicit
# binary directory, as the designing study's CMakeLists.txt does; the API's JSON is header-only.
#
#   cmake -S . -B ~/.cache/science-build/exp-district-economy-problem -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-district-economy-problem -j4
#   (from the repository's root)
#   ~/.cache/science-build/exp-district-economy-problem/calibrate research/experiments/2026-11-05-a-district-economy-the-problem/results/problem-inputs.json research/experiments/2026-11-05-a-district-economy-the-problem/results
#   ~/.cache/science-build/exp-district-economy-problem/check_economy research/experiments/2026-11-05-a-district-economy-the-problem/results/problem-inputs.json research/experiments/2026-11-05-a-district-economy-the-problem/results/problem.json

cmake_minimum_required(VERSION 3.20)
project(experiment_district_economy_problem LANGUAGES CXX)

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

# The header's own dependencies, for the studies to link: economy.hpp and the designing study's lanes.hpp and
# walking.hpp beside it.
add_library(district_economy INTERFACE)
target_include_directories(district_economy INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${SCIENCE}/api/cpp)
target_link_libraries(district_economy INTERFACE urban::events urban::landscape urban::affinity architecture::space
                      physis::statistical_mechanics physis::computational_geometry heuresis::variation calx::linalg calx::random)

add_executable(calibrate calibrate.cpp)
target_link_libraries(calibrate PRIVATE district_economy urban::rhythm urban::response Threads::Threads)
target_compile_options(calibrate PRIVATE -Wall -Wextra)

add_executable(check_economy check_economy.cpp)
target_link_libraries(check_economy PRIVATE district_economy Threads::Threads)
target_compile_options(check_economy PRIVATE -Wall -Wextra)
