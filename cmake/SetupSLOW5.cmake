include(${CMAKE_CURRENT_LIST_DIR}/Util.cmake)

function(add_slow5_to_target TARGET_NAME)
    if(NOSLOW5)
        target_compile_definitions(${TARGET_NAME} PRIVATE NSLOW5RH=1)
    else()
        if(SLOW5_COMPILE)
            add_dependencies(${TARGET_NAME} slow5_build)
        endif()
        link_imported_library(${TARGET_NAME} slow5 ${SLOW5_DIR})
    endif()
endfunction()

function(setup_slow5)
    if(NOT NOSLOW5)
        set(SLOW5_SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/slow5lib)
        if(SLOW5_COMPILE)
            if(NOT SLOW5_DIR)
                override_cached(SLOW5_DIR ${EXTERNAL_PROJECTS_BUILD_DIR}/slow5lib)
            endif()
            message(STATUS "Compiling slow5 to ${SLOW5_DIR}")
            ExternalProject_Add(
                slow5_build
                BINARY_DIR ${SLOW5_DIR}
                SOURCE_DIR ${SLOW5_SOURCE_DIR}
                INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory ${SLOW5_SOURCE_DIR}/include ${SLOW5_DIR}/include
                             && ${CMAKE_COMMAND} -E make_directory ${SLOW5_DIR}/lib
                             && ${CMAKE_COMMAND} -E rename ${SLOW5_DIR}/libslow5.so ${SLOW5_DIR}/lib/libslow5.so
            )
            message(STATUS "Current dir: ${CMAKE_CURRENT_BINARY_DIR}")
        else()
            if(NOT SLOW5_DIR)
                message(FATAL_ERROR "SLOW5_COMPILE is OFF, but no dir provided")
            endif()
        endif()
        message(STATUS "Using slow5 from ${SLOW5_DIR}")
        define_imported_library(slow5 ${SLOW5_DIR})
    endif()
endfunction()
