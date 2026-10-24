# mapping: a minimap2-like mapper. Minimizer anchors, chaining by dynamic programming with Li's
# gap cost, primary and secondary chains and their mapping quality, base-level alignment of the
# gaps between seeds and of the ends; SAM; evaluation against the simulator's truth.
add_library(genomics_mapping INTERFACE)
add_library(genomics::mapping ALIAS genomics_mapping)
target_include_directories(genomics_mapping INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_mapping INTERFACE cxx_std_20)
target_link_libraries(genomics_mapping INTERFACE genomics::index genomics::align genomics::sequence calx::scheduler)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_mapping_<topic>.
  file(GLOB mapping_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${mapping_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_mapping_${topic} SOURCES ${source} LIBRARIES genomics::mapping genomics::simulate calx::random)
  endforeach()
endif()
