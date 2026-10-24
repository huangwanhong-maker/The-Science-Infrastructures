# hmm: hidden Markov models (Viterbi, forward, backward, posteriors, Baum–Welch, sampling),
# the CpG-island model, and profile HMMs built from alignments.
add_library(genomics_hmm INTERFACE)
add_library(genomics::hmm ALIAS genomics_hmm)
target_include_directories(genomics_hmm INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(genomics_hmm INTERFACE genomics::sequence calx::random calx::scheduler)

if(GENOMICS_CHECKS)
  file(GLOB hmm_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${hmm_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_hmm_${topic} SOURCES ${source} LIBRARIES genomics::hmm genomics::simulate)
  endforeach()
endif()
