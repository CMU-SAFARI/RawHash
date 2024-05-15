include(${CMAKE_CURRENT_LIST_DIR}/Util.cmake)

function(setup_slow5)
    if(NOSLOW5)
        target_compile_definitions(${TARGET_NAME} PRIVATE NSLOW5RH=1)
    else()
        set(SLOW5_SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/slow5lib)
        if(SLOW5_COMPILE)
            if(NOT SLOW5_DIR)
                override_cached(SLOW5_DIR ${WORKDIR}/slow5lib)
            endif()
            ExternalProject_Add(
                slow5_build
                SOURCE_DIR ${SLOW5_SOURCE_DIR}
                BINARY_DIR ${SLOW5_DIR}/build
                INSTALL_COMMAND # slow5 doesn't have native install target
                    ${CMAKE_COMMAND} -E make_directory ${SLOW5_DIR}/lib &&
                    ${CMAKE_COMMAND} -E copy ${SLOW5_DIR}/build/libslow5.so ${SLOW5_DIR}/lib &&
                    ${CMAKE_COMMAND} -E copy_directory ${SLOW5_SOURCE_DIR}/include ${SLOW5_DIR}/include
            )
            add_dependencies(${TARGET_NAME} slow5_build)
        else()
            if(NOT SLOW5_DIR)
                message(FATAL_ERROR "SLOW5_COMPILE is OFF, but no dir provided")
            endif()
        endif()
        link_imported_library(slow5 ${SLOW5_DIR})
    endif()
endfunction()
