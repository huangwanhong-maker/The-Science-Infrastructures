# align: pairwise alignment. Biological scoring (substitution matrices, affine gaps, Karlin–Altschul
# statistics) on morphe's generic alignment over semirings, which defines the answers; the scalar
# schedules of the same sum (Gotoh, Myers–Miller, banded, X-drop), the fast kernels on vector lanes
# (bit-parallel, striped, wavefront, partial order), parametric alignment on morphe's polytope
# algebra, and the pair HMM.
add_library(genomics_align INTERFACE)
add_library(genomics::align ALIAS genomics_align)
target_include_directories(genomics_align INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(genomics_align INTERFACE genomics::sequence morphe::formal calx::cpu calx::solvers)

if(GENOMICS_CHECKS)
  # Every tests/check_<topic>.cpp is a check program, check_genomics_align_<topic>.
  file(GLOB align_checks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_*.cpp)
  foreach(source ${align_checks})
    get_filename_component(topic ${source} NAME_WE)
    string(REGEX REPLACE "^check_" "" topic ${topic})
    genomics_add_check(check_genomics_align_${topic} SOURCES ${source} LIBRARIES genomics::align genomics::simulate calx::random)
  endforeach()
endif()
