# Kepler's third law rediscovered from JPL's ephemeris: a standalone project on the foundations (calx, morphe,
# physis, tyche, heuresis) and the astronomy domain, added by relative path as heuresis adds its siblings.
#
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
#   build/kepler <path to de440s.bsp> results

cmake_minimum_required(VERSION 3.20)
project(keplers_third_law LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
find_package(Threads REQUIRED)

set(SCIENCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../..)

# science_add(<directory under the root> <a target it defines>): the library, added once.
function(science_add directory target)
  if(NOT TARGET ${target})
    get_filename_component(name ${directory} NAME)
    add_subdirectory(${SCIENCE_ROOT}/${directory} ${CMAKE_BINARY_DIR}/libraries/${name})
  endif()
endfunction()

science_add(foundations/calx calx::linalg)
science_add(foundations/morphe morphe::topology)
science_add(foundations/physis physis::units)
science_add(foundations/tyche tyche::probability)
science_add(foundations/heuresis heuresis::programs)
science_add(domains/astronomy astronomy::ephemerides)

add_executable(kepler kepler.cpp)
target_compile_options(kepler PRIVATE -Wall -Wextra)
target_link_libraries(kepler PRIVATE astronomy::ephemerides heuresis::programs heuresis::multiobjective calx::solvers Threads::Threads)
