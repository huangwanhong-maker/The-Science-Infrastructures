# The evolution domain's C++ (ADR 0043): phylogeny, population genetics and selection, on the
# genomics domain. Header-only; each module a directory with its own CMakeLists.txt and the
# target evolution::<module>, joining the build when it exists.
#
# The domain stands on genomics (sequences, alignments, variants) and uses calx (linear
# algebra, special functions, random streams, task graphs, the report of its checks) and
# morphe. Inside The Science Infrastructures they are its siblings, added here when no larger
# project has added them.

cmake_minimum_required(VERSION 3.20)
project(evolution LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(EVOLUTION_TOP_LEVEL ON)
else()
  set(EVOLUTION_TOP_LEVEL OFF)
endif()
option(EVOLUTION_CHECKS "Build the evolution domain's checks" ${EVOLUTION_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# The genomics domain beneath, which brings calx and morphe.
if(NOT TARGET genomics::sequence)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../genomics/CMakeLists.txt)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../genomics ${CMAKE_CURRENT_BINARY_DIR}/genomics)
  else()
    message(FATAL_ERROR "the evolution domain needs the genomics domain: build it inside The Science Infrastructures")
  endif()
endif()
if(EVOLUTION_CHECKS)
  enable_testing()
endif()

# evolution_add_check(<name> SOURCES <files> LIBRARIES <targets>): a check program on calx's
# verification, knowing where the domain's fixtures are.
set(EVOLUTION_DIR ${CMAKE_CURRENT_SOURCE_DIR})
function(evolution_add_check name)
  cmake_parse_arguments(ARG "" "" "SOURCES;LIBRARIES" ${ARGN})
  calx_add_check(${name} SOURCES ${ARG_SOURCES} LIBRARIES ${ARG_LIBRARIES})
  target_compile_definitions(${name} PRIVATE EVOLUTION_FIXTURES="${EVOLUTION_DIR}/tests/fixtures")
endfunction()

# The modules in their order (README §2).
foreach(module phylogeny population selection)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${module}/CMakeLists.txt)
    add_subdirectory(${module})
  endif()
endforeach()
