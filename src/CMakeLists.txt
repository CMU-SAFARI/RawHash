cmake_minimum_required(VERSION 3.10)
project(RawHash2Src LANGUAGES CXX)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    message(STATUS "ccache found: ${CCACHE_PROGRAM}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
else()
    message(STATUS "ccache not found, using default compiler")
endif()

find_program(MOLD_PATH mold)
if(MOLD_PATH)
    message(STATUS "mold found: ${MOLD_PATH}")
    set(CMAKE_EXE_LINKER_FLAGS "-B${MOLD_PATH}")
else()
    message(STATUS "mold not found, using default linker")
endif()

# TODO: Remove host-specific configs
site_name(HOSTNAME)
if(HOSTNAME STREQUAL "mikado1" OR HOSTNAME STREQUAL "mikado")
    set(CMAKE_DEBUG_TYPE "Debug")
    set(HDF5_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/extern/hdf5/build/include")
    set(HDF5_LIB_DIR "${CMAKE_SOURCE_DIR}/extern/hdf5/build/lib")
    set(POD5_DOWNLOAD OFF)
    set(NOPOD5 ON)
endif()
if(HOSTNAME STREQUAL "mikado")
    set(RUCLIENT_ENABLED ON)
    set(tsan "")
    set(asan "")
endif()

if(PYBINDING)
    message(STATUS "Building with Python binding support")
    set(TARGET_NAME rawhash_pybinding)
    # TODO: Non-existent in git?
    add_executable(${TARGET_NAME} rawhash_mapper.cpp)
else()
    set(TARGET_NAME rawhash2)
    add_executable(${TARGET_NAME} main.cpp)
endif()

set_target_properties(${TARGET_NAME} PROPERTIES CXX_STANDARD 11)

target_compile_options(${TARGET_NAME} PRIVATE
    -Wall
    -fopenmp
    -march=native
    -pthread
    -O3
)
target_compile_definitions(${TARGET_NAME} PRIVATE HAVE_KALLOC)

find_package(Threads REQUIRED)
target_link_libraries(${TARGET_NAME} PRIVATE Threads::Threads)

target_sources(${TARGET_NAME} PRIVATE
    bseq.c
    dtw.cpp
    hit.c
    kalloc.c
    kthread.c
    lchain.c
    revent.c
    rindex.c
    rmap.cpp
    roptions.c
    rseed.c
    rsig.c
    rsketch.c
    rutils.c
    sequence_until.c
)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -x c++")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(${TARGET_NAME} PRIVATE
        -g
        -fsanitize=address
    )
endif()

if(PROFILE)
    target_compile_options(${TARGET_NAME} PRIVATE
        -g
        -fno-omit-frame-pointer
    )
    target_compile_definitions(${TARGET_NAME} PRIVATE
        PROFILERH=1
    )
endif()

if(RUCLIENT_ENABLED)
    message(STATUS "ruclient enabled")
    set_target_properties(${TARGET_NAME} PROPERTIES CXX_STANDARD 20)
    target_compile_definitions(${TARGET_NAME} PRIVATE RUCLIENT_ENABLED)
    target_sources(${TARGET_NAME} PRIVATE rawhash_ruclient.cpp)
    # TODO: Place everything in subdirectories of CMAKE_SOURCE_DIR?
    target_link_directories(${TARGET_NAME} PRIVATE 
        {CMAKE_SOURCE_DIR}/../readuntil_fake/build/
        ~/.local/lib/
    )
    include_directories(
        {CMAKE_SOURCE_DIR}/../readuntil_fake/external/spdlog/include
        {CMAKE_SOURCE_DIR}/../readuntil_fake/include
    )

    # TODO: Any less painful way to do this?
    target_link_libraries(${TARGET_NAME} 
        ont_device_client_LIB
        grpc_utils_LIB
        utils_LIB
        ru_method_LIB
        data_client_LIB
        grpc_utils_LIB
        readuntil_grpc_proto
        grpc++_reflection
        grpc++
        grpc
        re2
        upb_json_lib
        upb_textformat_lib
        upb_collections_lib
        upb
        utf8_range_lib
        z
        cares
        gpr
        absl_random_distributions
        absl_random_seed_sequences
        absl_random_internal_pool_urbg
        absl_random_internal_randen
        absl_random_internal_randen_hwaes
        absl_random_internal_randen_hwaes_impl
        absl_random_internal_randen_slow
        absl_random_internal_platform
        absl_random_internal_seed_material
        absl_random_seed_gen_exception
        ssl
        crypto
        address_sorting
        protobuf
        absl_log_internal_check_op
        absl_leak_check
        absl_die_if_null
        absl_log_internal_conditions
        absl_log_internal_message
        absl_log_internal_nullguard
        absl_examine_stack
        absl_log_internal_format
        absl_log_internal_proto
        absl_log_internal_log_sink_set
        absl_log_sink
        absl_log_entry
        absl_flags
        absl_flags_internal
        absl_flags_marshalling
        absl_flags_reflection
        absl_flags_config
        absl_flags_program_name
        absl_flags_private_handle_accessor
        absl_flags_commandlineflag
        absl_flags_commandlineflag_internal
        absl_log_initialize
        absl_log_globals
        absl_log_internal_globals
        absl_raw_hash_set
        absl_hash
        absl_city
        absl_low_level_hash
        absl_hashtablez_sampler
        absl_statusor
        absl_status
        absl_cord
        absl_cordz_info
        absl_cord_internal
        absl_cordz_functions
        absl_exponential_biased
        absl_cordz_handle
        absl_crc_cord_state
        absl_crc32c
        absl_crc_internal
        absl_crc_cpu_detect
        absl_bad_optional_access
        absl_str_format_internal
        absl_strerror
        absl_synchronization
        absl_stacktrace
        absl_symbolize
        absl_debugging_internal
        absl_demangle_internal
        absl_graphcycles_internal
        absl_kernel_timeout_internal
        absl_malloc_internal
        absl_time
        absl_civil_time
        absl_time_zone
        absl_bad_variant_access
        absl_utf8_validity
        absl_strings
        absl_int128
        absl_string_view
        absl_throw_delegate
        absl_strings_internal
        absl_base
        absl_raw_logging_internal
        absl_log_severity
        absl_spinlock_wait
        dl
        m
        rt
        channel_data_LIB
        current_read_data_LIB
        channel_stats_LIB
        utils_LIB
        threadsafe_queue_LIB
        spdlogd
    )
else()
    message(STATUS "ruclient disabled")
endif()

set(WORKDIR ${CMAKE_SOURCE_DIR}/extern)

set(POD5_VERSION "0.2.2")
set(POD5_URLDIR "pod5-${POD5_VERSION}-${CMAKE_SYSTEM_NAME}")
set(POD5_REPO "https://github.com/nanoporetech/pod5-file-format")

if(DEFINED ENV{POD5_INCLUDE_DIR})
    set(POD5_DOWNLOAD 0)
else()
    set(POD5_INCLUDE_DIR "${WORKDIR}/${POD5_URLDIR}/include")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(POD5_LIB "lib64")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm)")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
            set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-gcc7-arm64.tar.gz")
            set(POD5_LIB "lib")
        else()
            set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-arm64.tar.gz")
        endif()
    else()
        set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-x64.tar.gz")
    endif()

    set(POD5_LIB_DIR "${WORKDIR}/${POD5_URLDIR}/${POD5_LIB}")
    set(POD5_LIBRARIES "${POD5_LIB_DIR}/libpod5_format.a"
                       "${POD5_LIB_DIR}/libarrow.a"
                       "${POD5_LIB_DIR}/libjemalloc_pic.a")
elseif(SYSTEM_NAME STREQUAL "Darwin")
    set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-osx-11.0-arm64.tar.gz")
    set(POD5_LIB_DIR "${WORKDIR}/${POD5_URLDIR}/lib")
    set(POD5_LIBRARIES "${POD5_LIB_DIR}/libpod5_format.a"
                       "${POD5_LIB_DIR}/libarrow.a")
elseif(SYSTEM_NAME STREQUAL "Windows_NT")
    set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-win")
endif()

set(HDF5_DIR "${WORKDIR}/hdf5/")
set(HDF5_BUILD_DIR "${HDF5_DIR}/build")
if(DEFINED ENV{HDF5_INCLUDE_DIR})
    set(HDF5_COMPILE 0)
else()
    set(HDF5_INCLUDE_DIR "${HDF5_BUILD_DIR}/include")
endif()
set(HDF5_LIB_DIR "${HDF5_BUILD_DIR}/lib")
set(HDF5_LIB "hdf5")

set(SLOW5_DIR "${WORKDIR}/slow5lib/")
if(DEFINED ENV{SLOW5_INCLUDE_DIR})
    set(SLOW5_COMPILE 0)
else()
    set(SLOW5_INCLUDE_DIR "${SLOW5_DIR}/include")
endif()
set(SLOW5_LIB_DIR "${SLOW5_DIR}/lib")

include_directories("${WORKDIR}/tensorflow")
target_link_directories(${TARGET_NAME} PRIVATE "${WORKDIR}/tensorflow/tflite_build/")

if(NOPOD5)
    target_compile_definitions(${TARGET_NAME} PRIVATE NPOD5RH=1)
else()
    include_directories("${POD5_INCLUDE_DIR}")
    target_link_directories(${TARGET_NAME} PRIVATE "${WORKDIR}/zstd/lib/")
    target_link_libraries(${TARGET_NAME} PRIVATE "${POD5_LIBRARIES}" "zstd")
endif()

if(NOHDF5)
    target_compile_definitions(${TARGET_NAME} PRIVATE NHDF5RH=1)
else()
    include_directories("${HDF5_INCLUDE_DIR}")
    target_link_libraries(${TARGET_NAME} PRIVATE "${HDF5_LIB_DIR}/lib${HDF5_LIB}.a")
endif()

if(NOSLOW5)
    target_compile_definitions(${TARGET_NAME} PRIVATE NSLOW5RH=1)
else()
    include_directories("${SLOW5_INCLUDE_DIR}")
    target_link_libraries(${TARGET_NAME} PRIVATE "${SLOW5_LIB_DIR}/libslow5.a")
endif()

target_link_libraries(${TARGET_NAME} PRIVATE "m" "z" "dl" "tensorflowlite_c")

if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
    target_compile_options(${TARGET_NAME} PRIVATE -D_FILE_OFFSET_BITS=64 -fsigned-char)
elseif(DEFINED ARM_NEON)
    target_compile_options(${TARGET_NAME} PRIVATE -D_FILE_OFFSET_BITS=64 -mfpu=neon -fsigned-char)
endif()

if(ENABLE_ASAN)
    message(STATUS "AddressSanitizer enabled")
    target_compile_options(${TARGET_NAME} PRIVATE -fsanitize=address)
    target_link_libraries(${TARGET_NAME} PRIVATE -fsanitize=address)
endif()

if(ENABLE_TSAN)
    message(STATUS "ThreadSanitizer enabled")
    target_compile_options(${TARGET_NAME} PRIVATE -fsanitize=thread)
    target_link_libraries(${TARGET_NAME} PRIVATE -fsanitize=thread)
endif()

if(NOT NOPOD5)
    add_custom_target(zstd
        COMMAND ${CMAKE_COMMAND} -E chdir ${WORKDIR}/zstd make -j
        COMMENT "Building ZSTD library"
    )
endif()

if(NOT NOHDF5 AND HDF5_COMPILE)
    add_custom_target(hdf5
        COMMAND ${CMAKE_COMMAND} -E echo "Configuring and building HDF5..."
        COMMAND ${CMAKE_COMMAND} -E make_directory ${HDF5_DIR}/build
        COMMAND ${CMAKE_COMMAND} -E chdir ${HDF5_DIR} ./configure --enable-threadsafe --disable-hl --prefix="${HDF5_BUILD_DIR}"
        COMMAND ${CMAKE_COMMAND} -E chdir ${HDF5_DIR} make -j
        COMMAND ${CMAKE_COMMAND} -E chdir ${HDF5_DIR} make install
    )
endif()

if(NOT NOSLOW5 AND SLOW5_COMPILE)
    add_custom_target(slow5
        COMMAND ${CMAKE_COMMAND} -E chdir ${SLOW5_DIR} make
    )
endif()

if(NOT NOPOD5 AND POD5_DOWNLOAD)
    message(STATUS "${CMAKE_COMMAND} -E chdir ${WORKDIR}/${POD5_URLDIR} wget -qO- ${POD5_URL} | tar -xzv")
    add_custom_target(pod5
        COMMAND ${CMAKE_COMMAND} -E make_directory ${WORKDIR}/${POD5_URLDIR}
        COMMAND wget -qO- ${POD5_URL} | tar -xzv -C ${WORKDIR}/${POD5_URLDIR}
        DEPENDS zstd
    )
endif()

add_custom_target(tensorflow
    COMMAND ${CMAKE_COMMAND} -E make_directory ${WORKDIR}/tensorflow/tflite_build
    COMMAND ${CMAKE_COMMAND} -E chdir ${WORKDIR}/tensorflow/tflite_build ${CMAKE_COMMAND} ../tensorflow/lite/c
    COMMAND ${CMAKE_COMMAND} -E chdir ${WORKDIR}/tensorflow/tflite_build ${CMAKE_COMMAND} --build . -j
)
