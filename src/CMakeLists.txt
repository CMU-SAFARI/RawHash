cmake_minimum_required(VERSION 3.10)
project(RawHash2Src)

option(USE_CCACHE "Use ccache to speed up rebuilds" ON)
option(USE_MOLD "Use mold linker for faster linking" ON)

option(PYBINDING "Build with Python bindings" OFF) # WIP
option(PROFILE "Enable profiling support" OFF)
option(RUCLIENT_ENABLED "Enable ReadUntil client" OFF)

option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(ENABLE_TSAN "Enable ThreadSanitizer" OFF)

option(POD5_DOWNLOAD "Download POD5" ON)
option(HDF5_COMPILE "Compile HDF5" ON)
option(SLOW5_COMPILE "Compile slow5lib" ON)

option(NOPOD5 "Do not use POD5" OFF)
option(NOHDF5 "Do not use HDF5" ON)
option(NOSLOW5 "Do not use SLOW5" OFF)

set(POD5_DIR "" CACHE PATH "Path to POD5 directory")
set(HDF5_DIR "" CACHE PATH "Path to HDF5 directory")
set(SLOW5_DIR "" CACHE PATH "Path to SLOW5 directory")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin CACHE PATH "Binaries directory")

set(WORKDIR ${CMAKE_BINARY_DIR}/extern)
file(MAKE_DIRECTORY ${WORKDIR})

set(TARGET_NAME rawhash2)

include(../cmake/SetupCCacheMold.cmake)
include(../cmake/SetupMainTarget.cmake)
include(../cmake/SetupRUClient.cmake)
include(../cmake/SetupPOD5.cmake)
include(../cmake/SetupHDF5.cmake)
include(../cmake/SetupSLOW5.cmake)
include(../cmake/SetupTFLite.cmake)

setup_ccache_mold()
setup_main_target(${TARGET_NAME})
setup_ruclient(${TARGET_NAME})
setup_pod5(${TARGET_NAME})
setup_hdf5(${TARGET_NAME})
setup_slow5(${TARGET_NAME})
setup_tflite(${TARGET_NAME})
