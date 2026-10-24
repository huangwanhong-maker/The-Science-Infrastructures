# sketch: MinHash bottom-s and FracMinHash sketches of k-mer sets, the Jaccard index, containment,
# the Mash distance and ANI; HyperLogLog; all-against-all distance matrices as task graphs.
add_library(genomics_sketch INTERFACE)
add_library(genomics::sketch ALIAS genomics_sketch)
target_include_directories(genomics_sketch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_sketch INTERFACE cxx_std_20)
target_link_libraries(genomics_sketch INTERFACE genomics::sequence calx::scheduler calx::special_functions)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_sketch_<topic>.
  file(GLOB sketch_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${sketch_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_sketch_${topic} SOURCES ${source} LIBRARIES genomics::sketch genomics::simulate)
  endforeach()
endif()
