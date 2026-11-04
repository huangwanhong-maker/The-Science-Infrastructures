# Fitting neurons, and the repertoire of their firing: a case study of heuresis on the neuroscience domain's
# spiking neurons (README.md). A project of its own, built against the foundations and the domain by relative
# path, each added once (as foundations/heuresis/CMakeLists.txt adds its siblings).
#
#   cmake -S . -B <build> -DCMAKE_BUILD_TYPE=Release
#   cmake --build <build> --target fitting_neurons -j4
#   <build>/fitting_neurons <this directory>/results

cmake_minimum_required(VERSION 3.20)
project(fitting_neurons LANGUAGES CXX)

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

# The foundations in their order, then the domain.
science_library(foundations/calx calx::linalg calx)
science_library(foundations/morphe morphe::topology morphe)
science_library(foundations/physis physis::units physis)
science_library(foundations/tyche tyche::probability tyche)
science_library(foundations/heuresis heuresis::search heuresis)
science_library(domains/neuroscience neuroscience::snn neuroscience)

add_executable(fitting_neurons fitting_neurons.cpp)
target_link_libraries(fitting_neurons PRIVATE heuresis::search heuresis::variation heuresis::strategies heuresis::swarms
                                              heuresis::diversity neuroscience::snn calx::scheduler Threads::Threads)
target_compile_options(fitting_neurons PRIVATE -Wall -Wextra -Wpedantic)
