# assembly: the bidirected de Bruijn graph of canonical k-mers and its unitigs, cleaning (tips,
# bubbles, low coverage) recorded step by step, contigs and their statistics, and the directed
# de Bruijn graph of a k-mer composition with Eulerian paths and their number by the BEST theorem.
add_library(genomics_assembly INTERFACE)
add_library(genomics::assembly ALIAS genomics_assembly)
target_include_directories(genomics_assembly INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(genomics_assembly INTERFACE genomics::kmers genomics::sequence morphe::combinatorics morphe::algebra morphe::arithmetic)

if(GENOMICS_CHECKS)
  file(GLOB assembly_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${assembly_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_assembly_${topic} SOURCES ${source} LIBRARIES genomics::assembly genomics::simulate)
  endforeach()
endif()
