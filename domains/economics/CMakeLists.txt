# The economics domain's C++: value carried between positions (transfer), its circulation and
# stocks (circulation), production as a circular flow (production), the spacetime complex of a
# record of transactions (spacetime), the faces a map gives a network of places (places), growth
# and the short run (growth, macro), and the partners and gravity of trade (trade).
# Header-only; each module a target economics::<module>.
#
# The domain uses morphe (complexes, the worldline construction), physis (cochains and the Hodge
# decomposition, the spherical Delaunay triangulation), calx (linear algebra, integrators and solvers,
# the report of its checks), tyche (regression) and nomos (measures of concentration and
# specialisation, Poisson regression with absorbed effects). Inside The
# Science Infrastructures they are the sibling foundations, added here when no larger project has
# added them.

cmake_minimum_required(VERSION 3.20)
project(economics LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(ECONOMICS_TOP_LEVEL ON)
else()
  set(ECONOMICS_TOP_LEVEL OFF)
endif()
option(ECONOMICS_CHECKS "Build the economics domain's checks" ${ECONOMICS_TOP_LEVEL})

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
    message(FATAL_ERROR "the economics domain needs calx, morphe and physis: build it inside Noetarium")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::discretisation)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(NOT TARGET tyche::regression)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/tyche ${CMAKE_CURRENT_BINARY_DIR}/tyche)
endif()
if(NOT TARGET nomos::measurement)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/nomos ${CMAKE_CURRENT_BINARY_DIR}/nomos)
endif()
if(ECONOMICS_CHECKS)
  enable_testing()
endif()

# transfer: transfer networks, their holonomies, loop deficits, circuits, prices and retention.
add_library(economics_transfer INTERFACE)
add_library(economics::transfer ALIAS economics_transfer)
target_include_directories(economics_transfer INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/transfer/include)
target_compile_features(economics_transfer INTERFACE cxx_std_20)
target_link_libraries(economics_transfer INTERFACE physis::discretisation morphe::topology morphe::combinatorics)

# circulation: the Hodge decomposition of flows, the pairing of flows and transports, balances by
# kind, and stocks under a stochastic rule.
add_library(economics_circulation INTERFACE)
add_library(economics::circulation ALIAS economics_circulation)
target_include_directories(economics_circulation INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/circulation/include)
target_link_libraries(economics_circulation INTERFACE economics::transfer calx::linalg)

# production: input-output tables, the Leontief model and its dual, the circulation of deliveries.
add_library(economics_production INTERFACE)
add_library(economics::production ALIAS economics_production)
target_include_directories(economics_production INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/production/include)
target_link_libraries(economics_production INTERFACE economics::circulation calx::linalg)

# spacetime: the spacetime complex of a record of transactions (morphe's worldline complex), its
# transports, curvature and reach.
add_library(economics_spacetime INTERFACE)
add_library(economics::spacetime ALIAS economics_spacetime)
target_include_directories(economics_spacetime INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/spacetime/include)
target_link_libraries(economics_spacetime INTERFACE economics::transfer morphe::topology)

# places: places on the sphere, and the faces their spherical Delaunay triangulation gives a network.
add_library(economics_places INTERFACE)
add_library(economics::places ALIAS economics_places)
target_include_directories(economics_places INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/places/include)
target_link_libraries(economics_places INTERFACE economics::transfer physis::computational_geometry)

# growth: the Solow–Swan model (Cobb–Douglas and CES), growth accounting and the perpetual inventory.
add_library(economics_growth INTERFACE)
add_library(economics::growth ALIAS economics_growth)
target_include_directories(economics_growth INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/growth/include)
target_link_libraries(economics_growth INTERFACE calx::integrators calx::solvers)

# macro: the IS–LM model of the short run, its multipliers and aggregate demand.
add_library(economics_macro INTERFACE)
add_library(economics::macro ALIAS economics_macro)
target_include_directories(economics_macro INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/macro/include)
target_compile_features(economics_macro INTERFACE cxx_std_20)

# trade: the partners of trade (concentration, intensity, two-way trade) and its gravity (PPML).
add_library(economics_trade INTERFACE)
add_library(economics::trade ALIAS economics_trade)
target_include_directories(economics_trade INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/trade/include)
target_link_libraries(economics_trade INTERFACE nomos::measurement nomos::econometrics tyche::regression physis::computational_geometry calx::linalg)

# regimes: the business cycle perceived from growth, a hidden Markov model fitted to many economies (ADR 0074).
add_library(economics_regimes INTERFACE)
add_library(economics::regimes ALIAS economics_regimes)
target_include_directories(economics_regimes INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/regimes/include)
target_link_libraries(economics_regimes INTERFACE tyche::processes)

# policy: decisions on the perceived regime, a policymaker's POMDP, a fund under a budget, a budget on support.
add_library(economics_policy INTERFACE)
add_library(economics::policy ALIAS economics_policy)
target_include_directories(economics_policy INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/policy/include)
target_link_libraries(economics_policy INTERFACE economics::regimes tyche::decision tyche::processes)

if(ECONOMICS_CHECKS)
  calx_add_check(check_economics_policy SOURCES policy/tests/check_policy.cpp LIBRARIES economics::policy economics::regimes calx::random)
  calx_add_check(check_economics_growth SOURCES growth/tests/check_growth.cpp LIBRARIES economics::growth economics::macro)
  calx_add_check(check_economics_places SOURCES places/tests/check_places.cpp LIBRARIES economics::places economics::circulation calx::random)
  calx_add_check(check_economics_transfer SOURCES transfer/tests/check_transfer.cpp LIBRARIES economics::transfer calx::random)
  calx_add_check(check_economics_circulation SOURCES circulation/tests/check_circulation.cpp LIBRARIES economics::circulation calx::random)
  calx_add_check(check_economics_production SOURCES production/tests/check_production.cpp LIBRARIES economics::production)
  calx_add_check(check_economics_trade SOURCES trade/tests/check_trade.cpp LIBRARIES economics::trade economics::circulation calx::random)
  calx_add_check(check_economics_gravity SOURCES trade/tests/check_gravity.cpp LIBRARIES economics::trade calx::random)
  calx_add_check(check_economics_spacetime SOURCES spacetime/tests/check_spacetime.cpp LIBRARIES economics::spacetime economics::circulation calx::random)
endif()
