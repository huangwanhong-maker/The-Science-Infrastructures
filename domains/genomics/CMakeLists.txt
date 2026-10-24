# The genomics domain's C++ (ADR 0043): sequences and their alphabets, simulation, k-mers,
# alignment, indexes, sketches, mapping, assembly, variants, hidden Markov models, multiple
# alignment, expression and annotation. Header-only; each module a directory with its own
# CMakeLists.txt and the target genomics::<module>, joining the build when it exists.
#
# The domain uses calx (lanes, task graphs, random streams, special functions, the report of
# its checks) and morphe (words and their indexes, alignment over semirings, exact integers).
# Inside The Science Infrastructures they are the sibling foundations, added here when no
# larger project has added them.

cmake_minimum_required(VERSION 3.20)
project(genomics LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(GENOMICS_TOP_LEVEL ON)
else()
  set(GENOMICS_TOP_LEVEL OFF)
endif()
option(GENOMICS_CHECKS "Build the genomics domain's checks" ${GENOMICS_TOP_LEVEL})
option(GENOMICS_NATIVE "Build the checks for this machine's widest lanes (-march=native)" ON)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(SCIENCE_FOUNDATIONS ${CMAKE_CURRENT_SOURCE_DIR}/../../foundations)
if(NOT TARGET calx::verification)
  if(EXISTS ${SCIENCE_FOUNDATIONS}/calx/CMakeLists.txt)
    add_subdirectory(${SCIENCE_FOUNDATIONS}/calx ${CMAKE_CURRENT_BINARY_DIR}/calx)
  else()
    message(FATAL_ERROR "the genomics domain needs calx and morphe: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::check)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(GENOMICS_CHECKS)
  enable_testing()
endif()

# genomics_add_check(<name> SOURCES <files> LIBRARIES <targets>): a check program on calx's
# verification, built for this machine's lanes where the compiler allows it, knowing where the
# domain's fixtures are and where the data library's cache is by default (SCIENCE_CACHE, when
# set, takes precedence at run time).
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-march=native GENOMICS_HAS_MARCH_NATIVE)
set(GENOMICS_DIR ${CMAKE_CURRENT_SOURCE_DIR})
function(genomics_add_check name)
  cmake_parse_arguments(ARG "" "" "SOURCES;LIBRARIES" ${ARGN})
  calx_add_check(${name} SOURCES ${ARG_SOURCES} LIBRARIES ${ARG_LIBRARIES})
  if(GENOMICS_NATIVE AND GENOMICS_HAS_MARCH_NATIVE)
    target_compile_options(${name} PRIVATE -march=native)
  endif()
  target_compile_definitions(${name} PRIVATE GENOMICS_FIXTURES="${GENOMICS_DIR}/tests/fixtures"
                                             GENOMICS_CACHE="${GENOMICS_DIR}/../../data/cache")
endfunction()

# The modules in their order (README §2).
foreach(module sequence simulate kmers align index sketch mapping assembly variant hmm msa expression annotation)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${module}/CMakeLists.txt)
    add_subdirectory(${module})
  endif()
endforeach()
