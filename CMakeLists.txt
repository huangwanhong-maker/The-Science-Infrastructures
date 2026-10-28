# The Science Infrastructures: the C++ of the foundations and the domains, in their order
# (docs/mathematics.md §10; ADR 0035).
#
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j
#   ctest --test-dir build --output-on-failure
#
# Each library is also a CMake project of its own; the runner of the checks
# (development/tools/run-tests.mjs) builds each one separately.

cmake_minimum_required(VERSION 3.20)
project(the_science_infrastructures LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(SCIENCE_CHECKS "Build the checks of every library" ON)
set(CALX_CHECKS ${SCIENCE_CHECKS})
set(MORPHE_CHECKS ${SCIENCE_CHECKS})
set(PHYSIS_CHECKS ${SCIENCE_CHECKS})
set(NEUROSCIENCE_CHECKS ${SCIENCE_CHECKS})
set(METROLOGY_CHECKS ${SCIENCE_CHECKS})
set(GEODESY_CHECKS ${SCIENCE_CHECKS})
set(ASTRONOMY_CHECKS ${SCIENCE_CHECKS})
set(GEOPHYSICS_CHECKS ${SCIENCE_CHECKS})
set(METEOROLOGY_CHECKS ${SCIENCE_CHECKS})
set(NAVIGATION_CHECKS ${SCIENCE_CHECKS})
set(HISTORY_CHECKS ${SCIENCE_CHECKS})
set(GENOMICS_CHECKS ${SCIENCE_CHECKS})
set(EVOLUTION_CHECKS ${SCIENCE_CHECKS})
set(ECONOMICS_CHECKS ${SCIENCE_CHECKS})
set(PLANETARY_SCIENCE_CHECKS ${SCIENCE_CHECKS})
if(SCIENCE_CHECKS)
  enable_testing()
endif()

# The foundations: calx and morphe stand alone; physis uses both.
add_subdirectory(foundations/calx)
add_subdirectory(foundations/morphe)
add_subdirectory(foundations/physis)

# The domains with C++, in the order of the domains.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/domains/neuroscience/CMakeLists.txt)
  add_subdirectory(domains/neuroscience)
endif()
add_subdirectory(domains/metrology)
add_subdirectory(domains/history)
add_subdirectory(domains/economics)
add_subdirectory(domains/genomics)
add_subdirectory(domains/geodesy)
add_subdirectory(domains/evolution)
add_subdirectory(domains/astronomy)
foreach(domain geophysics meteorology)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/domains/${domain}/CMakeLists.txt)
    add_subdirectory(domains/${domain})
  endif()
endforeach()
add_subdirectory(domains/navigation)
add_subdirectory(domains/planetary-science)
