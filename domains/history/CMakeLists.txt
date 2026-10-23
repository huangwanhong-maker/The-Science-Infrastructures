# The history domain's C++: the record and its chronology, the construction of the history
# complex, and its readings (ADR 0041). Header-only; each module a target history::<module>.
#
# The domain uses morphe (complexes, homology, π₁) and calx (the report of its checks). Inside
# The Science Infrastructures they are the sibling foundations, added here when no larger
# project has added them.

cmake_minimum_required(VERSION 3.20)
project(history LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(HISTORY_TOP_LEVEL ON)
else()
  set(HISTORY_TOP_LEVEL OFF)
endif()
option(HISTORY_CHECKS "Build the history domain's checks" ${HISTORY_TOP_LEVEL})

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
    message(FATAL_ERROR "the history domain needs calx and morphe: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_FOUNDATIONS}/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(HISTORY_CHECKS)
  enable_testing()
endif()

# record: the record of a history, read from JSON lines; its dates and calendars.
add_library(history_record INTERFACE)
add_library(history::record ALIAS history_record)
target_include_directories(history_record INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/record/include)
target_compile_features(history_record INTERFACE cxx_std_20)

# construction: the history complex of a record, and the loops of its H₁.
add_library(history_construction INTERFACE)
add_library(history::construction ALIAS history_construction)
target_include_directories(history_construction INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/construction/include)
target_link_libraries(history_construction INTERFACE history::record morphe::topology)

# readings: the signed ties as a ℤ₂ field in temporal gauge; the disciplines as colours.
add_library(history_readings INTERFACE)
add_library(history::readings ALIAS history_readings)
target_include_directories(history_readings INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/readings/include)
target_link_libraries(history_readings INTERFACE history::construction)

# foam: the refined complex, the spin foam of a history at G = ℤ₂, and the other complexes.
add_library(history_foam INTERFACE)
add_library(history::foam ALIAS history_foam)
target_include_directories(history_foam INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/foam/include)
target_link_libraries(history_foam INTERFACE history::readings calx::random)

if(HISTORY_CHECKS)
  # The canons are the domain's records; the checks read them by the path given here.
  set(HISTORY_CANONS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/canons)
  calx_add_check(check_history_chronology SOURCES record/tests/check_chronology.cpp LIBRARIES history::record)
  calx_add_check(check_history_record SOURCES record/tests/check_record.cpp LIBRARIES history::record)
  calx_add_check(check_history_construction SOURCES construction/tests/check_construction.cpp LIBRARIES history::construction)
  calx_add_check(check_history_readings SOURCES readings/tests/check_readings.cpp LIBRARIES history::readings)
  calx_add_check(check_history_foam SOURCES foam/tests/check_foam.cpp LIBRARIES history::foam)
  foreach(check check_history_record check_history_construction check_history_readings check_history_foam)
    target_compile_definitions(${check} PRIVATE HISTORY_CANONS="${HISTORY_CANONS_DIR}")
  endforeach()
endif()
