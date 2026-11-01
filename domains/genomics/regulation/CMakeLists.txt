# regulation: expression from regulatory sequence: one-hot, k-mer and gapped k-mer encodings; the
# sigma-70 promoter model learnt from transcription start sites and scanned along genomes;
# thermodynamic models of occupancy and fold change; and sequence-to-expression models learnt from
# measured activities (tyche's elastic net and convolutional networks), their filters read as motifs.
add_library(genomics_regulation INTERFACE)
add_library(genomics::regulation ALIAS genomics_regulation)
target_include_directories(genomics_regulation INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_regulation INTERFACE cxx_std_20)
target_link_libraries(genomics_regulation INTERFACE genomics::annotation genomics::sequence tyche::learning tyche::regression tyche::descriptive calx::random
                                                    calx::scheduler)

if(GENOMICS_CHECKS)
  file(GLOB regulation_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${regulation_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_regulation_${topic} SOURCES ${source} LIBRARIES genomics::regulation)
  endforeach()
endif()
