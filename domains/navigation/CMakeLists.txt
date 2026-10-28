# The navigation bridge's C++: the catalogue of frames of the Solar System (frames), joining
# astronomy's ephemerides and rotations, metrology's time scales and Earth orientation, and
# geodesy's figure of the Earth in one frame graph. Header-only; each module a target
# navigation::<module>.
#
# The bridge stands on astronomy, geodesy and geography (and through them metrology), and uses
# physis and calx. Its checks read the kernels and the Earth orientation from the infrastructure's
# data cache (jpl-de440s, naif-lunar-frames, iers-eop) and skip what is missing.

cmake_minimum_required(VERSION 3.20)
project(navigation LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(NAVIGATION_TOP_LEVEL ON)
else()
  set(NAVIGATION_TOP_LEVEL OFF)
endif()
option(NAVIGATION_CHECKS "Build the navigation bridge's checks" ${NAVIGATION_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(SCIENCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)
if(NOT TARGET calx::verification)
  if(EXISTS ${SCIENCE_ROOT}/foundations/calx/CMakeLists.txt)
    add_subdirectory(${SCIENCE_ROOT}/foundations/calx ${CMAKE_CURRENT_BINARY_DIR}/calx)
  else()
    message(FATAL_ERROR "the navigation bridge needs calx and physis: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_ROOT}/foundations/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::geometry)
  add_subdirectory(${SCIENCE_ROOT}/foundations/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
# The domains beneath, in their order, when no larger project has added them.
if(NOT TARGET metrology::frames)
  add_subdirectory(${SCIENCE_ROOT}/domains/metrology ${CMAKE_CURRENT_BINARY_DIR}/metrology)
endif()
if(NOT TARGET geodesy::figure)
  add_subdirectory(${SCIENCE_ROOT}/domains/geodesy ${CMAKE_CURRENT_BINARY_DIR}/geodesy)
endif()
if(NOT TARGET astronomy::rotation)
  add_subdirectory(${SCIENCE_ROOT}/domains/astronomy ${CMAKE_CURRENT_BINARY_DIR}/astronomy)
endif()
if(NAVIGATION_CHECKS)
  enable_testing()
endif()

# orbits: two-line element sets, the SGP4/SDP4 propagator and the TEME frame.
add_library(navigation_orbits INTERFACE)
add_library(navigation::orbits ALIAS navigation_orbits)
target_include_directories(navigation_orbits INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/orbits/include)
target_compile_features(navigation_orbits INTERFACE cxx_std_20)
target_link_libraries(navigation_orbits INTERFACE metrology::frames metrology::time physis::geometry)

# frames: the catalogue of the Solar System's frames (and the Galaxy's), places on the Earth and the Moon,
# and addresses: the language of points of spacetime, read and answered in any frame and chart.
add_library(navigation_frames INTERFACE)
add_library(navigation::frames ALIAS navigation_frames)
target_include_directories(navigation_frames INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/frames/include)
target_compile_features(navigation_frames INTERFACE cxx_std_20)
target_link_libraries(navigation_frames INTERFACE metrology::frames geodesy::figure geodesy::projection astronomy::rotation physis::mechanics navigation::orbits)

if(NAVIGATION_CHECKS)
  calx_add_check(check_navigation_sgp4 SOURCES orbits/tests/check_sgp4.cpp LIBRARIES navigation::orbits)
  target_compile_definitions(check_navigation_sgp4 PRIVATE NAVIGATION_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures")
  calx_add_check(check_navigation_catalogue SOURCES frames/tests/check_catalogue.cpp LIBRARIES navigation::frames)
  target_compile_definitions(check_navigation_catalogue PRIVATE NAVIGATION_ROOT="${CMAKE_CURRENT_SOURCE_DIR}")
  calx_add_check(check_navigation_address SOURCES frames/tests/check_address.cpp LIBRARIES navigation::frames)
  target_compile_definitions(check_navigation_address PRIVATE NAVIGATION_ROOT="${CMAKE_CURRENT_SOURCE_DIR}"
                             NAVIGATION_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures")
endif()

# vehicle: a vehicle's frames: the aerospace Euler angles, body rates, ground speed, track and flight
# path, air data, the body frame on the Earth.
add_library(navigation_vehicle INTERFACE)
add_library(navigation::vehicle ALIAS navigation_vehicle)
target_include_directories(navigation_vehicle INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/vehicle/include)
target_compile_features(navigation_vehicle INTERFACE cxx_std_20)
target_link_libraries(navigation_vehicle INTERFACE geodesy::figure physis::geometry)

if(NAVIGATION_CHECKS)
  calx_add_check(check_navigation_vehicle SOURCES vehicle/tests/check_attitude.cpp LIBRARIES navigation::vehicle calx::random)
endif()
