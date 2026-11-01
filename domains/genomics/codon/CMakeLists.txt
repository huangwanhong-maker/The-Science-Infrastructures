# codon: codon usage and expression: the coding and tRNA genes of an annotated genome, synonymous
# families under any genetic code, effective numbers of codons, CAI, Fop, MILC, Karlin and Mrázek's
# predicted expression, the tRNA adaptation index with its S value and stAI, and the agreement of
# an index with measured protein abundances.
add_library(genomics_codon INTERFACE)
add_library(genomics::codon ALIAS genomics_codon)
target_include_directories(genomics_codon INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_codon INTERFACE cxx_std_20)
target_link_libraries(genomics_codon INTERFACE genomics::annotation genomics::sequence tyche::descriptive calx::optimisation calx::random calx::scheduler
                                               calx::special_functions)

if(GENOMICS_CHECKS)
  file(GLOB codon_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${codon_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_codon_${topic} SOURCES ${source} LIBRARIES genomics::codon)
  endforeach()
endif()
