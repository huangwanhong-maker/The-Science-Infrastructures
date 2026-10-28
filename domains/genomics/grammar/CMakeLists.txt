# grammar: group-equivariant grammars of nucleic acids (the user's GE-CSG draft, made exact). The
# Klein group of the bases, its subgroups (IUPAC's codes) and characters (the classifications of
# the bases), breakings; the strand group on words, canonical k-mers by Burnside's lemma, and the
# symmetries a genome's statistics keep; hairpin grammars and their symmetry; the genetic code as
# a function on (C₂)⁶. The general machinery is morphe's (algebra/action, formal/equivariant).
add_library(genomics_grammar INTERFACE)
add_library(genomics::grammar ALIAS genomics_grammar)
target_include_directories(genomics_grammar INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(genomics_grammar INTERFACE genomics::sequence morphe::formal morphe::algebra)

if(GENOMICS_CHECKS)
  file(GLOB grammar_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${grammar_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_grammar_${topic} SOURCES ${source} LIBRARIES genomics::grammar calx::random)
  endforeach()
endif()
