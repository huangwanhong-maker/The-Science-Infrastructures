# A walk past Melbourne's sensors: the travelling-salesman problem on the City of Melbourne's pedestrian
# counters, solved by the search foundation heuresis on geodesic distances from the geodesy domain.
#
# A standalone project: the foundations and the domain are added from the infrastructure by relative path,
# each once (as foundations/heuresis/CMakeLists.txt adds its siblings).

cmake_minimum_required(VERSION 3.20)
project(walk_past_melbournes_sensors LANGUAGES CXX)

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
if(NOT TARGET geodesy::geodesic)
  add_subdirectory(${SCIENCE}/domains/geodesy ${CMAKE_BINARY_DIR}/geodesy)
endif()

add_executable(walk walk.cpp)
target_compile_options(walk PRIVATE -Wall -Wextra)
target_link_libraries(walk PRIVATE heuresis::genetic heuresis::swarms heuresis::local heuresis::annealing heuresis::distributions heuresis::analysis
                                   geodesy::geodesic geodesy::figure)
