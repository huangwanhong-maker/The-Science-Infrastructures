# karyotype: the description of chromosomes, ISCN (karyotypes, arr, seq), parsed by its grammar on
# morphe's BNF reader; the human ideogram's bands as a tree; karyotypes interpreted as derivative
# chromosomes over oriented segments, copy number and fusions; arr and seq descriptions written
# from copy-number segments and adjacencies.
add_library(genomics_karyotype INTERFACE)
add_library(genomics::karyotype ALIAS genomics_karyotype)
target_include_directories(genomics_karyotype INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_karyotype INTERFACE cxx_std_20)
target_link_libraries(genomics_karyotype INTERFACE genomics::sequence morphe::formal)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_karyotype_<topic>.
  file(GLOB karyotype_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${karyotype_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_karyotype_${topic} SOURCES ${source} LIBRARIES genomics::karyotype)
  endforeach()
endif()
