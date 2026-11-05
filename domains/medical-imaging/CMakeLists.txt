# The medical-imaging domain's C++: images of the body and the structures reconstructed from them (ADR 0070).
# Header-only; each module a target medical_imaging::<module>:
#   dicom           DICOM files read (uncompressed, little endian), series sorted and checked, their grid in the
#                   patient frame LPS; NIfTI's RAS carried to it;
#   modalities      CT's Hounsfield units, windows and forward model; TOF-MRA's inflow signal and forward model;
#   reconstruction  vessels from one image: detected by vesselness, delineated at the local half-maximum,
#                   thinned to a vascular tree with radii, bounded by a surface;
#   integration     vessels from several images: registered by mutual information, fused, and told apart as
#                   arteries and veins by what each modality shows.
#
# The domain stands on anatomy (frames, vessels) and the foundations calx, morphe, physis and eikon
# (module.json); eikon holds every method on images, the domain only what is medical.

cmake_minimum_required(VERSION 3.20)
project(medical_imaging LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(MEDICAL_IMAGING_TOP_LEVEL ON)
else()
  set(MEDICAL_IMAGING_TOP_LEVEL OFF)
endif()
option(MEDICAL_IMAGING_CHECKS "Build the medical-imaging domain's checks" ${MEDICAL_IMAGING_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(SCIENCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)
if(NOT TARGET calx::verification)
  if(EXISTS ${SCIENCE_ROOT}/foundations/calx/CMakeLists.txt)
    add_subdirectory(${SCIENCE_ROOT}/foundations/calx ${CMAKE_CURRENT_BINARY_DIR}/calx)
  else()
    message(FATAL_ERROR "the medical-imaging domain needs the foundations: build it inside Noetarium")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_ROOT}/foundations/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::computational_geometry)
  add_subdirectory(${SCIENCE_ROOT}/foundations/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(NOT TARGET eikon::grid)
  add_subdirectory(${SCIENCE_ROOT}/foundations/eikon ${CMAKE_CURRENT_BINARY_DIR}/eikon)
endif()
if(NOT TARGET anatomy::vessels)
  add_subdirectory(${SCIENCE_ROOT}/domains/anatomy ${CMAKE_CURRENT_BINARY_DIR}/anatomy)
endif()
if(MEDICAL_IMAGING_CHECKS)
  enable_testing()
endif()

function(medical_imaging_module name)
  cmake_parse_arguments(ARG "" "" "LIBRARIES;CHECK_LIBRARIES" ${ARGN})
  add_library(medical_imaging_${name} INTERFACE)
  add_library(medical_imaging::${name} ALIAS medical_imaging_${name})
  target_include_directories(medical_imaging_${name} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${name}/include)
  target_compile_features(medical_imaging_${name} INTERFACE cxx_std_20)
  target_link_libraries(medical_imaging_${name} INTERFACE ${ARG_LIBRARIES})
  if(MEDICAL_IMAGING_CHECKS)
    calx_add_check(check_medical_imaging_${name} SOURCES ${name}/tests/check_${name}.cpp LIBRARIES medical_imaging::${name} calx::random ${ARG_CHECK_LIBRARIES})
  endif()
endfunction()

medical_imaging_module(dicom LIBRARIES anatomy::frames eikon::grid)
medical_imaging_module(modalities LIBRARIES eikon::filters eikon::phantoms calx::random)
medical_imaging_module(reconstruction LIBRARIES anatomy::vessels eikon::features eikon::segmentation eikon::skeleton eikon::surfaces
                       CHECK_LIBRARIES medical_imaging::modalities eikon::fusion eikon::phantoms)
medical_imaging_module(integration LIBRARIES medical_imaging::reconstruction eikon::fusion eikon::registration
                       CHECK_LIBRARIES medical_imaging::modalities eikon::phantoms)

# The reconstruction on real CT (the Medical Segmentation Decathlon's hepatic vessels, from the data library's
# product; skipped when it is not built).
if(MEDICAL_IMAGING_CHECKS)
  calx_add_check(check_medical_imaging_msd SOURCES reconstruction/tests/check_msd.cpp LIBRARIES medical_imaging::reconstruction eikon::fusion)
  target_compile_definitions(check_medical_imaging_msd PRIVATE MEDICAL_IMAGING_PRODUCTS="${SCIENCE_ROOT}/data/cache/products")
endif()
