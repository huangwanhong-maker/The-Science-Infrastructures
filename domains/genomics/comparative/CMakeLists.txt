# comparative: whole-genome alignment of two genomes by maximal unique matches (MUMmer's method):
# the MUMs of both strands from morphe's suffix array and LCP array, their heaviest chain by
# morphe's increasing subsequences, the gaps between anchors, collinear blocks, the dot plot.
add_library(genomics_comparative INTERFACE)
add_library(genomics::comparative ALIAS genomics_comparative)
target_include_directories(genomics_comparative INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_comparative INTERFACE cxx_std_20)
target_link_libraries(genomics_comparative INTERFACE genomics::sequence morphe::formal morphe::combinatorics)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_comparative_<topic>.
  file(GLOB comparative_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${comparative_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_comparative_${topic} SOURCES ${source} LIBRARIES genomics::comparative genomics::simulate calx::random)
  endforeach()
endif()
