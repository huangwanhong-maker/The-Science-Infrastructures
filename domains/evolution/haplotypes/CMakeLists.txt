# haplotypes: Li and Stephens's copying model of a haplotype on a reference panel (forward,
# backward, posteriors and Viterbi in O(N K)), the pair model of a diploid's genotypes (its
# posterior, phasing as the joint maximum, switch errors), imputation from both, and the
# haplotype scans for selection (EHH, iHH, iHS standardised within bins of frequency).
add_library(evolution_haplotypes INTERFACE)
add_library(evolution::haplotypes ALIAS evolution_haplotypes)
target_include_directories(evolution_haplotypes INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(evolution_haplotypes INTERFACE evolution::population genomics::hmm calx::random calx::scheduler)

if(EVOLUTION_CHECKS)
  file(GLOB haplotypes_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${haplotypes_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    evolution_add_check(check_evolution_haplotypes_${topic} SOURCES ${source} LIBRARIES evolution::haplotypes)
  endforeach()
endif()
