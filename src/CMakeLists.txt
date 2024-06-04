cmake_minimum_required(VERSION 3.10)
project(RawHash2Src)

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
set(RUCLIENT_DIR "" CACHE PATH "Path to ReadUntil directory (already built or where it should be built)")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin CACHE PATH "Directory where to put binaries")

set(EXTERNAL_PROJECTS_BUILD_DIR ${CMAKE_BINARY_DIR}/extern) # default build directory for external projects
file(MAKE_DIRECTORY ${EXTERNAL_PROJECTS_BUILD_DIR})
message(STATUS "External projects build directory: ${EXTERNAL_PROJECTS_BUILD_DIR}")

include(../cmake/SetupCCacheMold.cmake)
include(../cmake/SetupRawHashLikeTarget.cmake)
include(../cmake/SetupRUClient.cmake)
# include(../cmake/SetupPOD5.cmake)
include(../cmake/SetupHDF5.cmake)
include(../cmake/SetupSLOW5.cmake)
include(../cmake/SetupTFLite.cmake)

# setup_slow5(${TARGET_NAME})

set(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "verbose makefile" FORCE) # print compilation commands

enable_ccache()
set_alternate_linker(mold)

# todo: currently can only compile rawhash2 or the below, not both
# this happens because functions like setup_hdf5 define the hdf5 target and the relationship to the target name
# rather than doing it separately

set(TARGET_NAME rawhash2)
add_executable(${TARGET_NAME} main.cpp)
setup_rawhashlike_target(${TARGET_NAME})
# setup_pod5(${TARGET_NAME}) # does not seem to work when pod5 is disabled
target_compile_definitions(${TARGET_NAME} PRIVATE NPOD5RH=1) # todo: workaround
setup_ruclient(${TARGET_NAME})
setup_hdf5(${TARGET_NAME})
setup_slow5(${TARGET_NAME})
setup_tflite(${TARGET_NAME})

# # set(TARGET_NAME rawhash2_wrapper_example)
# # add_executable(${TARGET_NAME} wrapper_example.cpp rawhash_wrapper.cpp rawhash_wrapper.hpp)
# set(TARGET_NAME rawhash2_wrapper)
# add_library(${TARGET_NAME} SHARED rawhash_wrapper.cpp rawhash_wrapper.hpp)
# setup_rawhashlike_target(${TARGET_NAME})
# # setup_pod5(${TARGET_NAME})
# target_compile_definitions(${TARGET_NAME} PRIVATE NPOD5RH=1) # todo: workaround
# setup_ruclient(${TARGET_NAME})
# setup_hdf5(${TARGET_NAME})
# setup_slow5(${TARGET_NAME})
# setup_tflite(${TARGET_NAME})

# set(TARGET_NAME rawhash2_usinglib)
# add_executable(${TARGET_NAME} wrapper_example.cpp)
# target_link_libraries(${TARGET_NAME} PRIVATE rawhash2_wrapper)
# # setup_rawhashlike_target(${TARGET_NAME})
# # # setup_pod5(${TARGET_NAME})
# # target_compile_definitions(${TARGET_NAME} PRIVATE NPOD5RH=1) # todo: workaround
# # setup_ruclient(${TARGET_NAME})
# # setup_hdf5(${TARGET_NAME})
# # setup_slow5(${TARGET_NAME})
# # setup_tflite(${TARGET_NAME})


# setup_pod5()
# add_pod5_to_target(${TARGET_NAME})
