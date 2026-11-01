# somatic: the somatic genome of a tumour against its normal (ADR 0064): MuTect's tumour and normal
# tests and filters, logR and BAF with ASPCF and ASCAT's purity, ploidy and allele-specific copy
# number, cancer cell fractions and their clones, catalogues and signatures, MSIsensor's test.
add_library(genomics_somatic INTERFACE)
add_library(genomics::somatic ALIAS genomics_somatic)
target_include_directories(genomics_somatic INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_somatic INTERFACE cxx_std_20)
target_link_libraries(genomics_somatic INTERFACE genomics::variant genomics::mapping genomics::sequence tyche::time_series tyche::regression tyche::multivariate
                                                 calx::special_functions calx::scheduler calx::solvers calx::random calx::linalg)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_somatic_<topic>.
  file(GLOB somatic_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${somatic_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_somatic_${topic} SOURCES ${source} LIBRARIES genomics::somatic genomics::simulate calx::random)
  endforeach()
endif()
