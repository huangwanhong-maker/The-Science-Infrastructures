# expression: transcriptomes and simulated RNA-seq, pseudoalignment by k-mer equivalence classes,
# abundances by EM, and differential expression by negative binomial models with Wald tests and
# the false discovery rate.
add_library(genomics_expression INTERFACE)
add_library(genomics::expression ALIAS genomics_expression)
target_include_directories(genomics_expression INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(genomics_expression INTERFACE genomics::sequence genomics::kmers calx::special_functions calx::linalg calx::solvers calx::random
                                                    calx::scheduler)

if(GENOMICS_CHECKS)
  file(GLOB expression_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${expression_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_expression_${topic} SOURCES ${source} LIBRARIES genomics::expression genomics::simulate)
  endforeach()
endif()
