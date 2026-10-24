# population: genotype matrices and allele counts, Hardy–Weinberg tests, linkage disequilibrium,
# diversity, spectra and F_ST, the coalescent with growth and recombination, Wright–Fisher
# populations forward in time, and population structure (principal components, admixture).
add_library(evolution_population INTERFACE)
add_library(evolution::population ALIAS evolution_population)
target_include_directories(evolution_population INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(evolution_population INTERFACE calx::random calx::linalg calx::special_functions calx::quadrature
                                                      calx::scheduler calx::cpu)
target_link_libraries(evolution_population INTERFACE genomics::variant)

if(EVOLUTION_CHECKS)
  # The checks are built for this machine (popcount for the packed haplotypes, the widest lanes
  # for the kinship matrix), as the genomics domain's are.
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-march=native EVOLUTION_POPULATION_MARCH_NATIVE)
  file(GLOB population_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${population_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    evolution_add_check(check_evolution_population_${topic} SOURCES ${source} LIBRARIES evolution::population)
    if(EVOLUTION_POPULATION_MARCH_NATIVE)
      target_compile_options(check_evolution_population_${topic} PRIVATE -march=native)
    endif()
  endforeach()
endif()
