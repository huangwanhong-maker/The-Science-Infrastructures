# annotation: features of GenBank and GFF3 as tracks; position weight matrices, their
# information and score distributions, Gibbs sampling and MEME's EM; a self-trained gene
# finder; exact repeats by the suffix array.
add_library(genomics_annotation INTERFACE)
add_library(genomics::annotation ALIAS genomics_annotation)
target_include_directories(genomics_annotation INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(genomics_annotation INTERFACE genomics::sequence genomics::hmm morphe::formal calx::random calx::scheduler calx::special_functions)

if(GENOMICS_CHECKS)
  file(GLOB annotation_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${annotation_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_annotation_${topic} SOURCES ${source} LIBRARIES genomics::annotation genomics::simulate)
  endforeach()
endif()
