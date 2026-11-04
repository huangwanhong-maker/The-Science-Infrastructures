# A district's economy as a field theory: the people of the district as a field on its space–time complex, whose law
# is the stationary point and the Euclidean measure of an action S; the plan as a field coupled to it; and searches
# that use the action's structure, against black-box searches at equal budgets. The third of three studies of one
# shared problem (../2026-11-05-a-district-economy-the-problem).
#
# A standalone project: the foundations (calx, morphe, physis, tyche, heuresis) and the domains (architecture, then
# urban science) are added from the infrastructure by relative path, each once with an explicit binary directory,
# as the shared problem's CMakeLists.txt does; the problem's directory is added for its INTERFACE target
# district_economy (its own programs are left out of the build).
#
#   cmake -S . -B ~/.cache/science-build/exp-economy-field -DCMAKE_BUILD_TYPE=Release
#   cmake --build ~/.cache/science-build/exp-economy-field -j4
#   (from the repository's root)
#   ~/.cache/science-build/exp-economy-field/check_field <problem-inputs.json> <problem.json> <sensors.json>
#   ~/.cache/science-build/exp-economy-field/field_theory <problem-inputs.json> <problem.json> <sensors.json> <counters.json> <results directory>

cmake_minimum_required(VERSION 3.20)
project(experiment_district_economy_field_theory LANGUAGES CXX)

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
if(NOT TARGET district_economy)
  # The shared problem: its INTERFACE target district_economy (economy.hpp); its calibrate and check_economy are not built.
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../2026-11-05-a-district-economy-the-problem ${CMAKE_BINARY_DIR}/the-problem EXCLUDE_FROM_ALL)
endif()

set(FIELD_LIBRARIES district_economy physis::discretisation physis::field_theory heuresis::search heuresis::genetic heuresis::annealing
                    heuresis::distributions heuresis::strategies heuresis::variation calx::transforms calx::linalg calx::random calx::scheduler
                    morphe::topology Threads::Threads)

add_executable(field_theory field_theory.cpp)
target_include_directories(field_theory PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(field_theory PRIVATE ${FIELD_LIBRARIES})
target_compile_options(field_theory PRIVATE -Wall -Wextra)

add_executable(check_field check_field.cpp)
target_include_directories(check_field PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(check_field PRIVATE ${FIELD_LIBRARIES})
target_compile_options(check_field PRIVATE -Wall -Wextra)
