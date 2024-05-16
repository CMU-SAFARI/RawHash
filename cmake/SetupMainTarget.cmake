function(setup_main_target TARGET_NAME)
    if(PYBINDING)
        message(FATAL_ERROR "Building with Python binding support is not implemented")
    else()
        add_executable(${TARGET_NAME} main.cpp)
    endif()
    set_target_properties(${TARGET_NAME} PROPERTIES CXX_STANDARD 11)
    set_target_properties(${TARGET_NAME} PROPERTIES C_STANDARD 11)

    find_package(Threads REQUIRED)
    target_link_libraries(${TARGET_NAME} PRIVATE Threads::Threads m z dl)
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
        kalloc.c
        kthread.c
        revent.c
        rmap.cpp
        roptions.c
        rsketch.c
        rutils.c
        sequence_until.c
    )
    # C files that rely on hdf5_tools.hpp
    # Should be compiled as CXX for now
    set(PSEUDO_C_SOURCES
        hit.c
        lchain.c
        rindex.c
        rseed.c
        rsig.c
    )
    foreach(source IN LISTS PSEUDO_C_SOURCES)
        set_source_files_properties(${source} PROPERTIES LANGUAGE CXX)
    endforeach()
    target_sources(${TARGET_NAME} PRIVATE ${PSEUDO_C_SOURCES})


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
