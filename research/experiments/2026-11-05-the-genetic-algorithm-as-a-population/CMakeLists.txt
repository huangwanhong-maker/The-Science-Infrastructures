# The genetic algorithm as a population: a case study of heuresis beside the evolution domain's population
# genetics (README.md). A project of its own, built against the foundations and the domains by relative path,
# each added once (as foundations/heuresis/CMakeLists.txt adds its siblings).
#
#   cmake -S . -B <build> -DCMAKE_BUILD_TYPE=Release
#   cmake --build <build> --target ga_population -j4
#   <build>/ga_population <this directory>/results

cmake_minimum_required(VERSION 3.20)
project(ga_population LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
find_package(Threads REQUIRED)

set(SCIENCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../..)

# science_library(<path from the root> <a target it defines> <binary directory>): the library, added once.
function(science_library path target binary)
  if(NOT TARGET ${target})
    add_subdirectory(${SCIENCE_ROOT}/${path} ${CMAKE_BINARY_DIR}/${binary})
  endif()
endfunction()

# The foundations in their order, then the domains (evolution stands on genomics).
science_library(foundations/calx calx::linalg calx)
science_library(foundations/morphe morphe::topology morphe)
science_library(foundations/physis physis::units physis)
science_library(foundations/tyche tyche::probability tyche)
science_library(foundations/heuresis heuresis::search heuresis)
science_library(domains/genomics genomics::sequence genomics)
science_library(domains/evolution evolution::population evolution)

add_executable(ga_population ga_population.cpp)
target_link_libraries(ga_population PRIVATE heuresis::search heuresis::variation heuresis::selection heuresis::genetic heuresis::benchmarks
                                            evolution::population calx::linalg calx::special_functions calx::scheduler Threads::Threads)
target_compile_options(ga_population PRIVATE -Wall -Wextra -Wpedantic)
