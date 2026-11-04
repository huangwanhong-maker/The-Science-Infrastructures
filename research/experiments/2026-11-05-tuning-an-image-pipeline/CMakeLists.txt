# Tuning an image pipeline: a case study of heuresis on eikon's reconstruction of tubes from images of
# synthetic phantoms (README.md). A project of its own, built against the foundations by relative path, each
# added once (as foundations/heuresis/CMakeLists.txt adds its siblings).
#
#   cmake -S . -B <build> -DCMAKE_BUILD_TYPE=Release
#   cmake --build <build> --target tune_pipeline -j4
#   <build>/tune_pipeline <this directory>/results

cmake_minimum_required(VERSION 3.20)
project(tune_pipeline LANGUAGES CXX)

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

# The foundations in their order.
science_library(foundations/calx calx::linalg calx)
science_library(foundations/morphe morphe::topology morphe)
science_library(foundations/physis physis::units physis)
science_library(foundations/tyche tyche::probability tyche)
science_library(foundations/heuresis heuresis::search heuresis)
science_library(foundations/eikon eikon::grid eikon)

add_executable(tune_pipeline tune_pipeline.cpp)
target_link_libraries(tune_pipeline PRIVATE heuresis::search heuresis::variation heuresis::strategies heuresis::multiobjective
                                            eikon::phantoms eikon::filters eikon::features eikon::segmentation eikon::skeleton eikon::fusion
                                            calx::scheduler Threads::Threads)
target_compile_options(tune_pipeline PRIVATE -Wall -Wextra -Wpedantic)
