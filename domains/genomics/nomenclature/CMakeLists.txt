# nomenclature: the description of sequence variants. HGVS (g., m., o., c., n., r., p.) parsed by its
# grammar on morphe's BNF reader and written back; changes described by the 3' rule with VCF's and
# SPDI's canonical forms; transcripts from GenBank and GFF3, the projection g. -> c./n. -> p. and the
# consequences of the Sequence Ontology; VCF's structural-variant and breakend records (ADR 0064).
add_library(genomics_nomenclature INTERFACE)
add_library(genomics::nomenclature ALIAS genomics_nomenclature)
target_include_directories(genomics_nomenclature INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_nomenclature INTERFACE cxx_std_20)
target_link_libraries(genomics_nomenclature INTERFACE genomics::sequence genomics::annotation genomics::variant morphe::formal)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_nomenclature_<topic>.
  file(GLOB nomenclature_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${nomenclature_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_nomenclature_${topic} SOURCES ${source} LIBRARIES genomics::nomenclature calx::random)
  endforeach()
endif()
