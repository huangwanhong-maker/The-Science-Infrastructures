# classify: taxonomic classification of reads by exact k-mers and the lowest common ancestor
# (Kraken's method): the taxonomy of the lineages of GenBank records, the database of canonical
# k-mers and their LCAs, the heaviest root-to-leaf path of a read's hits.
add_library(genomics_classify INTERFACE)
add_library(genomics::classify ALIAS genomics_classify)
target_include_directories(genomics_classify INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_classify INTERFACE cxx_std_20)
target_link_libraries(genomics_classify INTERFACE genomics::kmers genomics::sequence calx::scheduler)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_classify_<topic>.
  file(GLOB classify_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${classify_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_classify_${topic} SOURCES ${source} LIBRARIES genomics::classify genomics::simulate calx::random)
  endforeach()
endif()
