# structural: structural variants and copy number from paired reads (ADR 0064). Breakends and
# adjacencies; the library and the classes of pairs; discordant pairs, split reads, gaps and clips
# as evidence; LUMPY's breakpoints; events and SVTyper's genotypes; read depth with GC and
# mappability, its segmentation by tyche's change points, copy number and aneuploidy; B-allele
# fractions and bcftools roh's runs of homozygosity; evaluation as Truvari's; the calls described in
# VCF, HGVS (genomics::nomenclature) and ISCN (genomics::karyotype).
add_library(genomics_structural INTERFACE)
add_library(genomics::structural ALIAS genomics_structural)
target_include_directories(genomics_structural INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_structural INTERFACE cxx_std_20)
target_link_libraries(genomics_structural INTERFACE genomics::variant genomics::mapping genomics::hmm genomics::align genomics::sequence genomics::nomenclature genomics::karyotype tyche::time_series calx::scheduler)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_structural_<topic>.
  file(GLOB structural_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${structural_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_structural_${topic} SOURCES ${source} LIBRARIES genomics::structural genomics::simulate calx::random)
  endforeach()
endif()
