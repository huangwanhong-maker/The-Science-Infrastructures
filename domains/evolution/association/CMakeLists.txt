# association: genotype to trait and to expression — cohorts with truth, association scans, LD and
# clumping, cis-eQTL mapping with permutations, fine-mapping (ABF, SuSiE), colocalisation, expression
# predicted from genotypes (PrediXcan, S-PrediXcan), SNP heritability (GREML, Haseman–Elston) and
# polygenic scores (C+T, LDpred-inf). Its statistics are tyche's (regression, testing).
add_library(evolution_association INTERFACE)
add_library(evolution::association ALIAS evolution_association)
target_include_directories(evolution_association INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(evolution_association INTERFACE evolution::population evolution::haplotypes tyche::regression tyche::testing
                                                       calx::random calx::linalg calx::scheduler calx::cpu)

if(EVOLUTION_CHECKS)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-march=native EVOLUTION_ASSOCIATION_MARCH_NATIVE)
  file(GLOB association_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${association_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    evolution_add_check(check_evolution_association_${topic} SOURCES ${source} LIBRARIES evolution::association)
    if(EVOLUTION_ASSOCIATION_MARCH_NATIVE)
      target_compile_options(check_evolution_association_${topic} PRIVATE -march=native)
    endif()
  endforeach()
endif()
