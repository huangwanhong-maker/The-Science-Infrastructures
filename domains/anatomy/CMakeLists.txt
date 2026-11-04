# The anatomy domain's C++: the structures of the body (ADR 0066).
# Header-only; each module a target anatomy::<module>:
#   frames     the anatomical directions, the patient frames LPS and RAS, the cardinal planes and the
#              radiological and neurological displays of a section;
#   parts      FMA's concepts, the is-a and part-of relations as forests and their closures, the elements that
#              represent the concepts, the systems of the body by an ordered table, sides;
#   body       BodyParts3D's adult man read from the data library's product: the atlas and the meshes;
#   measures   a part's topology (morphe), its area, solid, principal axes and extents (physis), compounds,
#              and the parts that touch;
#   sections   the body cut by a cardinal plane in the display's coordinates; vertebral levels.
#
# The domain stands on the foundations alone (module.json): calx (verification), morphe (topology) and
# physis (computational_geometry, mechanics). The checks that read the body need the data library's product
# (data/cache/products/body-bp3d-*), and skip when it has not been built.

cmake_minimum_required(VERSION 3.20)
project(anatomy LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(ANATOMY_TOP_LEVEL ON)
else()
  set(ANATOMY_TOP_LEVEL OFF)
endif()
option(ANATOMY_CHECKS "Build the anatomy domain's checks" ${ANATOMY_TOP_LEVEL})

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
    message(FATAL_ERROR "the anatomy domain needs the foundations: build it inside The Science Infrastructures")
  endif()
endif()
if(NOT TARGET morphe::topology)
  add_subdirectory(${SCIENCE_ROOT}/foundations/morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
endif()
if(NOT TARGET physis::computational_geometry)
  add_subdirectory(${SCIENCE_ROOT}/foundations/physis ${CMAKE_CURRENT_BINARY_DIR}/physis)
endif()
if(ANATOMY_CHECKS)
  enable_testing()
endif()

function(anatomy_module name)
  cmake_parse_arguments(ARG "" "" "LIBRARIES" ${ARGN})
  add_library(anatomy_${name} INTERFACE)
  add_library(anatomy::${name} ALIAS anatomy_${name})
  target_include_directories(anatomy_${name} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${name}/include)
  target_compile_features(anatomy_${name} INTERFACE cxx_std_20)
  target_link_libraries(anatomy_${name} INTERFACE ${ARG_LIBRARIES})
  if(ANATOMY_CHECKS)
    calx_add_check(check_anatomy_${name} SOURCES ${name}/tests/check_${name}.cpp LIBRARIES anatomy::${name})
    target_compile_definitions(check_anatomy_${name} PRIVATE ANATOMY_PRODUCTS="${SCIENCE_ROOT}/data/cache/products"
                               ANATOMY_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures")
  endif()
endfunction()

anatomy_module(frames)
anatomy_module(parts)
anatomy_module(body LIBRARIES anatomy::parts physis::computational_geometry)
anatomy_module(measures LIBRARIES anatomy::body morphe::topology physis::computational_geometry physis::mechanics)
anatomy_module(sections LIBRARIES anatomy::frames anatomy::measures physis::computational_geometry)
