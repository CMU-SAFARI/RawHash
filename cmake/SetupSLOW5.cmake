include(ExternalProject)

function(setup_slow5)
    if(NOSLOW5)
        target_compile_definitions(${TARGET_NAME} PRIVATE NSLOW5RH=1)
    else()
        set(SLOW5_DIR ${WORKDIR}/slow5lib)
        set(SLOW5_BUILD_DIR ${SLOW5_DIR}/build)
        set(SLOW5_SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/slow5lib)
        if(SLOW5_COMPILE)
            if(NOT SLOW5_INCLUDE_DIR)
                override_cached(SLOW5_INCLUDE_DIR ${SLOW5_SOURCE_DIR}/include)
            endif()
            ExternalProject_Add(
                slow5_build
                SOURCE_DIR ${SLOW5_SOURCE_DIR}
                BINARY_DIR ${SLOW5_BUILD_DIR}
                INSTALL_DIR ${SLOW5_DIR}
                INSTALL_COMMAND ""
            )
            add_dependencies(${TARGET_NAME} slow5_build)
        else()
            if(NOT SLOW5_INCLUDE_DIR)
                message(FATAL_ERROR "SLOW5_COMPILE is OFF, but no include dir provided")
            endif()
        endif()
        add_library(slow5 STATIC IMPORTED)
        set_target_properties(slow5 PROPERTIES
            IMPORTED_LOCATION $${SLOW5_DIR}/lib/libslow5.a
            INTERFACE_INCLUDE_DIRECTORIES ${SLOW5_INCLUDE_DIR}
        )
        target_link_libraries(${TARGET_NAME} PRIVATE slow5)
    endif()
endfunction()
