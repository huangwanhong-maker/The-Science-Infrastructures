# simulate: genomes with repeats, diploid mutation with its truth, short and long reads, and
# families of sequences with their true alignment, all from seeded streams.
add_library(genomics_simulate INTERFACE)
add_library(genomics::simulate ALIAS genomics_simulate)
target_include_directories(genomics_simulate INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(genomics_simulate INTERFACE genomics::sequence calx::random)

if(GENOMICS_CHECKS)
  genomics_add_check(check_genomics_simulate SOURCES tests/check_simulate.cpp LIBRARIES genomics::simulate)
endif()
