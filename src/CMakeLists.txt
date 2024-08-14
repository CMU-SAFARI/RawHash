cmake_minimum_required(VERSION 3.10)
project(RawHash2)

# option(USE_CCACHE "Use ccache to speed up rebuilds" ON)
# option(USE_MOLD "Use mold linker for faster linking" ON)

# option(PYBINDING "Build with Python bindings" OFF) # WIP

# for debugging
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(ENABLE_TSAN "Enable ThreadSanitizer" OFF)
option(PROFILE "Enable profiling support" OFF)

# whether to compile or instead assume it is prebuilt
option(POD5_DOWNLOAD "Whether to build POD5" ON) # todo: rename to POD5_COMPILE
option(HDF5_COMPILE "Compile HDF5" ON)
option(SLOW5_COMPILE "Compile slow5lib" ON)

# features to enable/disable in RawHash
# option(NOPOD5 "Disable POD5 in build" OFF) # todo: figure out how to enable shared libs
option(NOPOD5 "Disable POD5 in build" ON)
option(NOHDF5 "Disable HDF5 in build" OFF)
option(NOSLOW5 "Disable SLOW5 in build" OFF)
option(RUCLIENT_ENABLED "Enable ReadUntil client" OFF)

set(POD5_DIR "" CACHE PATH "Path to POD5 directory (already built or where it should be built)")
set(HDF5_DIR "" CACHE PATH "Path to HDF5 directory (already built or where it should be built)")
set(SLOW5_DIR "" CACHE PATH "Path to SLOW5 directory (already built or where it should be built)")
set(RUCLIENT_DIR "" CACHE PATH "Path to ReadUntil directory (where it should be built)")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin CACHE PATH "Directory where to put binaries")

set(EXTERNAL_PROJECTS_BUILD_DIR ${CMAKE_BINARY_DIR}/extern) # default build directory for external projects
file(MAKE_DIRECTORY ${EXTERNAL_PROJECTS_BUILD_DIR})
message(STATUS "External projects build directory: ${EXTERNAL_PROJECTS_BUILD_DIR}")

include(${CMAKE_SOURCE_DIR}/cmake/SetupCCacheMold.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/SetupRawHashLikeTarget.cmake)

set(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "verbose makefile" FORCE) # print compilation commands

enable_ccache()
set_alternate_linker(mold)

add_executable(rawhash2_builtin main.cpp)
setup_rawhashlike_target(rawhash2_builtin)

add_library(rawhash2_lib SHARED rawhash_wrapper.cpp rawhash_wrapper.hpp)
setup_rawhashlike_target(rawhash2_lib)

file(GLOB_RECURSE HEADER_FILES ${CMAKE_SOURCE_DIR}/src/*.h ${CMAKE_SOURCE_DIR}/src/*.hpp)
set_property(TARGET rawhash2_lib PROPERTY PUBLIC_HEADER ${HEADER_FILES})
set_target_properties(rawhash2_lib PROPERTIES INSTALL_RPATH $ORIGIN OUTPUT_NAME rawhash2)
install(TARGETS rawhash2_lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/${PROJECT_NAME})

add_executable(rawhash2 wrapper_example.cpp)
target_link_libraries(rawhash2 PRIVATE rawhash2_lib)
set_target_properties(rawhash2 PROPERTIES INSTALL_RPATH $ORIGIN/../lib)
install(TARGETS rawhash2 RUNTIME DESTINATION bin)
