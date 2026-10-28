# signal: the spectra of DNA. The indicator sequences of the four bases and their transforms by
# calx's FFT, the spectrum and its mean, and the period-3 signal of protein-coding DNA, exactly in
# integers, over a sequence or in windows along a genome.
add_library(genomics_signal INTERFACE)
add_library(genomics::signal ALIAS genomics_signal)
target_include_directories(genomics_signal INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_signal INTERFACE cxx_std_20)
target_link_libraries(genomics_signal INTERFACE genomics::sequence calx::transforms)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_signal_<topic>.
  file(GLOB signal_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${signal_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_signal_${topic} SOURCES ${source} LIBRARIES genomics::signal genomics::simulate calx::random)
  endforeach()
endif()
