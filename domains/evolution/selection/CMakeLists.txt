# selection: synonymous and non-synonymous sites and differences by Nei and Gojobori, dN/dS
# under any NCBI genetic code; the McDonald–Kreitman test with Fisher's exact and the G-test.
add_library(evolution_selection INTERFACE)
add_library(evolution::selection ALIAS evolution_selection)
target_include_directories(evolution_selection INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(evolution_selection INTERFACE genomics::sequence calx::special_functions)

if(EVOLUTION_CHECKS)
  file(GLOB selection_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${selection_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    evolution_add_check(check_evolution_selection_${topic} SOURCES ${source} LIBRARIES evolution::selection)
  endforeach()
endif()
