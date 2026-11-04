# A district's economy as a flow: study (2) of the shared problem (../2026-11-05-a-district-economy-the-problem). People
# as a conserved flow on the district's streets, with its circulation; the search uses the model's structure.
#
# A standalone project: the foundations (calx, morphe, physis, tyche, nomos, heuresis) and the domains (architecture,
# urban science, economics) are added from the infrastructure by relative path, each once with an explicit binary
# directory; the shared problem's directory is added for its INTERFACE target `district_economy` (its own programs are
# left out of the build).
#
#   cmake -S . -B ~/.cache/science-build/exp-economy-flow -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-economy-flow -j4
#   (from the repository's root)
#   ~/.cache/science-build/exp-economy-flow/check_flow   <problem-inputs.json> <problem.json> <sensors.json>
#   ~/.cache/science-build/exp-economy-flow/as_a_flow    <problem-inputs.json> <problem.json> <sensors.json> <results directory>

cmake_minimum_required(VERSION 3.20)
project(experiment_district_economy_as_a_flow LANGUAGES CXX)

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
if(NOT TARGET nomos::measurement)
  add_subdirectory(${SCIENCE}/foundations/nomos ${CMAKE_BINARY_DIR}/nomos)
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
if(NOT TARGET economics::circulation)
  add_subdirectory(${SCIENCE}/domains/economics ${CMAKE_BINARY_DIR}/economics)
endif()
if(NOT TARGET district_economy)
  add_subdirectory(${PROBLEM} ${CMAKE_BINARY_DIR}/problem EXCLUDE_FROM_ALL)
endif()

set(FLOW_LIBRARIES district_economy urban::flow urban::rhythm economics::circulation physis::discretisation physis::computational_geometry tyche::regression
                   heuresis::genetic heuresis::strategies heuresis::search heuresis::variation calx::linalg calx::random calx::scheduler Threads::Threads)

add_executable(as_a_flow as_a_flow.cpp)
target_include_directories(as_a_flow PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(as_a_flow PRIVATE ${FLOW_LIBRARIES})
target_compile_options(as_a_flow PRIVATE -Wall -Wextra)

add_executable(check_flow check_flow.cpp)
target_include_directories(check_flow PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(check_flow PRIVATE ${FLOW_LIBRARIES})
target_compile_options(check_flow PRIVATE -Wall -Wextra)
