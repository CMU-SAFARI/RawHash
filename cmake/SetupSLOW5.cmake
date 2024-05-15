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
                BUILD_DIR ${SLOW5_DIR}
                SOURCE_DIR ${SLOW5_SOURCE_DIR}
                CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory ${SLOW5_SOURCE_DIR} ${SLOW5_DIR}
                BUILD_COMMAND make -C ${SLOW5_DIR}
                INSTALL_COMMAND ""
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
