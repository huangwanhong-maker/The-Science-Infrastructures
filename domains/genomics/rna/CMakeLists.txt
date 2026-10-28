# rna: RNA secondary structure by grammars. Pair tables and dot-bracket; Nussinov's maximal base
# pairing and the maximum expected accuracy structure by one recursion; Knudsen and Hein's
# stochastic context-free grammar and its reversal-closed presentation, as morphe's grammars,
# folded by morphe's best parse and inside–outside; EM training, untied or tied along the rule
# orbits of a group (genomics::grammar's symmetries of the bases and strands), and the
# likelihood-ratio test of a tie.
add_library(genomics_rna INTERFACE)
add_library(genomics::rna ALIAS genomics_rna)
target_include_directories(genomics_rna INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(genomics_rna INTERFACE genomics::sequence genomics::grammar morphe::formal morphe::algebra calx::random calx::scheduler
                                             calx::special_functions)

if(GENOMICS_CHECKS)
  file(GLOB rna_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${rna_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_rna_${topic} SOURCES ${source} LIBRARIES genomics::rna)
  endforeach()
endif()
