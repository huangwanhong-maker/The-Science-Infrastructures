# phylogeny: trees and Newick, splits and consensus; distances, UPGMA, neighbour joining and
# BIONJ; substitution models with Γ rates; Felsenstein's likelihood, branch lengths, NNI search
# and model parameters; parsimony; the bootstrap on task graphs; simulation along trees; and
# root-to-tip dating.
add_library(evolution_phylogeny INTERFACE)
add_library(evolution::phylogeny ALIAS evolution_phylogeny)
target_include_directories(evolution_phylogeny INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(evolution_phylogeny INTERFACE genomics::sequence calx::linalg calx::optimisation calx::special_functions calx::random
                                                     calx::scheduler calx::cpu)

if(EVOLUTION_CHECKS)
  # Built for this machine's widest lanes, which the pruning loops use.
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-march=native EVOLUTION_PHYLOGENY_MARCH_NATIVE)
  file(GLOB phylogeny_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${phylogeny_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    evolution_add_check(check_evolution_phylogeny_${topic} SOURCES ${source} LIBRARIES evolution::phylogeny)
    if(EVOLUTION_PHYLOGENY_MARCH_NATIVE)
      target_compile_options(check_evolution_phylogeny_${topic} PRIVATE -march=native)
    endif()
  endforeach()
endif()
