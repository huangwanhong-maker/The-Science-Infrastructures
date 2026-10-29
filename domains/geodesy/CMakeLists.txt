# The geodesy domain's C++: the figure of the Earth (ellipsoids, geodetic and Earth-centred
# coordinates, local level frames), geodesics (Karney 2013), map projections and grids (the
# transverse Mercator, UTM, UPS, MGRS), the geoid (EGM96's grid, heights above mean sea level) and normal
# gravity (the level ellipsoid of WGS 84 and GRS 80, its field at any height, and that field as one of
# physis's fields). Header-only; each module a target geodesy::<module>. The geoid's checks read the data
# library's grid from the cache (products/geoid/egm96-15.pgm) and skip it when it is missing.
#
# The domain uses physis (rotations and rigid motions, physis::geometry; fields, physis::fields) and calx
# (quadrature, and the report of its checks); it may use metrology, beneath it in the order of the
# domains. Inside The Science Infrastructures they are the sibling foundations and domains, added here
# when no larger project has added them.

cmake_minimum_required(VERSION 3.20)
project(geodesy LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(GEODESY_TOP_LEVEL ON)
else()
  set(GEODESY_TOP_LEVEL OFF)
endif()
option(GEODESY_CHECKS "Build the geodesy domain's checks" ${GEODESY_TOP_LEVEL})

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
    message(FATAL_ERROR "the geodesy domain needs calx and physis: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::geometry)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(GEODESY_CHECKS)
  enable_testing()
endif()

# figure: reference ellipsoids, geodetic and Earth-centred coordinates, local level frames.
add_library(geodesy_figure INTERFACE)
add_library(geodesy::figure ALIAS geodesy_figure)
target_include_directories(geodesy_figure INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/figure/include)
target_compile_features(geodesy_figure INTERFACE cxx_std_20)
target_link_libraries(geodesy_figure INTERFACE physis::geometry)

if(GEODESY_CHECKS)
  calx_add_check(check_geodesy_figure SOURCES figure/tests/check_figure.cpp LIBRARIES geodesy::figure)
  target_compile_definitions(check_geodesy_figure PRIVATE GEODESY_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures")
endif()

# geodesic: geodesics on the ellipsoid (Karney 2013): the direct and inverse problems, points along a line.
add_library(geodesy_geodesic INTERFACE)
add_library(geodesy::geodesic ALIAS geodesy_geodesic)
target_include_directories(geodesy_geodesic INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/geodesic/include)
target_compile_features(geodesy_geodesic INTERFACE cxx_std_20)
target_link_libraries(geodesy_geodesic INTERFACE geodesy::figure)

# projection: the transverse Mercator by Krüger's series (Karney 2011), UTM and UPS, MGRS.
add_library(geodesy_projection INTERFACE)
add_library(geodesy::projection ALIAS geodesy_projection)
target_include_directories(geodesy_projection INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/projection/include)
target_compile_features(geodesy_projection INTERFACE cxx_std_20)
target_link_libraries(geodesy_projection INTERFACE geodesy::figure)

# geoid: a grid of undulations (EGM96 at 15'), interpolated; heights above mean sea level.
add_library(geodesy_geoid INTERFACE)
add_library(geodesy::geoid ALIAS geodesy_geoid)
target_include_directories(geodesy_geoid INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/geoid/include)
target_compile_features(geodesy_geoid INTERFACE cxx_std_20)

# gravity: the level ellipsoid (WGS 84, GRS 80), its derived constants and normal gravity at any height,
# and normal gravity as a field (field.hpp).
add_library(geodesy_gravity INTERFACE)
add_library(geodesy::gravity ALIAS geodesy_gravity)
target_include_directories(geodesy_gravity INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/gravity/include)
target_compile_features(geodesy_gravity INTERFACE cxx_std_20)
target_link_libraries(geodesy_gravity INTERFACE geodesy::figure physis::fields calx::quadrature)

if(GEODESY_CHECKS)
  calx_add_check(check_geodesy_gravity SOURCES gravity/tests/check_gravity.cpp LIBRARIES geodesy::gravity calx::random)
  target_compile_definitions(check_geodesy_gravity PRIVATE GEODESY_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures")
  calx_add_check(check_geodesy_geodesic SOURCES geodesic/tests/check_geodesic.cpp LIBRARIES geodesy::geodesic)
  target_compile_definitions(check_geodesy_geodesic PRIVATE GEODESY_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures")
  calx_add_check(check_geodesy_projection SOURCES projection/tests/check_projection.cpp LIBRARIES geodesy::projection)
  target_compile_definitions(check_geodesy_projection PRIVATE GEODESY_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures")
  calx_add_check(check_geodesy_geoid SOURCES geoid/tests/check_geoid.cpp LIBRARIES geodesy::geoid)
  target_compile_definitions(check_geodesy_geoid PRIVATE GEODESY_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures"
                             GEODESY_GEOID_GRID="${CMAKE_CURRENT_SOURCE_DIR}/../../data/cache/products/geoid/egm96-15.pgm")
endif()
