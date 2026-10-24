# variant: pileups of mapped reads, genotype likelihoods and the Bayesian caller of a diploid,
# the normalisation of indels, VCF 4.2, and the evaluation of calls against a truth set.
add_library(genomics_variant INTERFACE)
add_library(genomics::variant ALIAS genomics_variant)
target_include_directories(genomics_variant INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_variant INTERFACE cxx_std_20)
target_link_libraries(genomics_variant INTERFACE genomics::mapping genomics::align genomics::sequence calx::special_functions calx::scheduler)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_variant_<topic>.
  file(GLOB variant_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${variant_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_variant_${topic} SOURCES ${source} LIBRARIES genomics::variant genomics::simulate calx::random)
  endforeach()
endif()
