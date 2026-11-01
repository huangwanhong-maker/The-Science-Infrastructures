# simulate: genomes with repeats, diploid mutation with its truth, short and long reads, and
# families of sequences with their true alignment, all from seeded streams.
add_library(genomics_simulate INTERFACE)
add_library(genomics::simulate ALIAS genomics_simulate)
target_include_directories(genomics_simulate INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(genomics_simulate INTERFACE genomics::sequence calx::random)

if(GENOMICS_CHECKS)
  genomics_add_check(check_genomics_simulate SOURCES tests/check_simulate.cpp LIBRARIES genomics::simulate)
  # Every other tests/check_<topic>.cpp is a check program, check_genomics_simulate_<topic>.
  file(GLOB simulate_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  list(REMOVE_ITEM simulate_checks ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_simulate.cpp)
  foreach(source ${simulate_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_simulate_${topic} SOURCES ${source} LIBRARIES genomics::simulate)
  endforeach()
endif()
