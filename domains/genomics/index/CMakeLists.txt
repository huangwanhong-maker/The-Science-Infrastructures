# index: the FMD-index of a genome and its reverse complement on morphe's SA-IS, with
# bidirectional extension and super-maximal exact matches; the minimizer index of a genome.
add_library(genomics_index INTERFACE)
add_library(genomics::index ALIAS genomics_index)
target_include_directories(genomics_index INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_index INTERFACE cxx_std_20)
target_link_libraries(genomics_index INTERFACE genomics::sequence morphe::formal calx::scheduler calx::random)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_index_<topic>.
  file(GLOB index_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${index_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_index_${topic} SOURCES ${source} LIBRARIES genomics::index genomics::simulate)
  endforeach()
endif()
