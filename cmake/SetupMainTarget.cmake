function(setup_main_target)
    if(PYBINDING)
        message(STATUS "Building with Python binding support")
        set(TARGET_NAME rawhash_pybinding CACHE INTERNAL "Main target")
        # TODO: Non-existent in git?
        add_executable(${TARGET_NAME} rawhash_mapper.cpp)
    else()
        set(TARGET_NAME rawhash2 CACHE INTERNAL "Main target")
        add_executable(${TARGET_NAME} main.cpp)
    endif()
    set_target_properties(${TARGET_NAME} PROPERTIES CXX_STANDARD 11)
    set_target_properties(${TARGET_NAME} PROPERTIES C_STANDARD 11)

    find_package(Threads REQUIRED)
    target_link_libraries(${TARGET_NAME} PRIVATE Threads::Threads)
    target_compile_options(${TARGET_NAME} PRIVATE -Wall -fopenmp -march=native -O3)
    target_compile_definitions(${TARGET_NAME} PRIVATE HAVE_KALLOC)

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
endfunction()
