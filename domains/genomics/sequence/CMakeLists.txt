# sequence: alphabets, packed DNA, genetic codes, FASTA, FASTQ and GenBank, k-mers and their
# samples, composition.
add_library(genomics_sequence INTERFACE)
add_library(genomics::sequence ALIAS genomics_sequence)
target_include_directories(genomics_sequence INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(genomics_sequence INTERFACE cxx_std_20)

if(GENOMICS_CHECKS)
  genomics_add_check(check_genomics_sequence SOURCES tests/check_sequence.cpp LIBRARIES genomics::sequence calx::random)
  genomics_add_check(check_genomics_formats SOURCES tests/check_formats.cpp LIBRARIES genomics::sequence calx::random)
  genomics_add_check(check_genomics_kmers SOURCES tests/check_kmers.cpp LIBRARIES genomics::sequence calx::random)
endif()
