# msa: multiple alignment. The alignment and its formats, guide trees (k-mer distances, UPGMA,
# CLUSTAL W weights), progressive alignment of profiles with MUSCLE's refinement, the centre
# star, alignment by partial orders, and the sum-of-pairs, Q and TC scores.
add_library(genomics_msa INTERFACE)
add_library(genomics::msa ALIAS genomics_msa)
target_include_directories(genomics_msa INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(genomics_msa INTERFACE genomics::sequence genomics::align calx::scheduler)

if(GENOMICS_CHECKS)
  file(GLOB msa_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${msa_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_msa_${topic} SOURCES ${source} LIBRARIES genomics::msa genomics::simulate calx::random)
  endforeach()
endif()
