# kmers: canonical k-mers counted in parallel and deterministically, the k-mer spectrum and what
# it says of the genome (error threshold, coverage, size, heterozygosity), Bloom filters and
# count-min sketches.
add_library(genomics_kmers INTERFACE)
add_library(genomics::kmers ALIAS genomics_kmers)
target_include_directories(genomics_kmers INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(genomics_kmers INTERFACE genomics::sequence calx::scheduler calx::random)

if(GENOMICS_CHECKS)
  file(GLOB kmers_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${kmers_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_kmers_${topic} SOURCES ${source} LIBRARIES genomics::kmers genomics::simulate)
  endforeach()
endif()
