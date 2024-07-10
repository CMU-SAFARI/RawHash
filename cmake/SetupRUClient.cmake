include(${CMAKE_CURRENT_LIST_DIR}/Util.cmake)

function(add_ruclient_to_target TARGET_NAME)
    if(RUCLIENT_ENABLED)
        set_target_properties(${TARGET_NAME} PROPERTIES CXX_STANDARD 20)
        target_compile_definitions(${TARGET_NAME} PRIVATE RUCLIENT_ENABLED)
        target_sources(${TARGET_NAME} PRIVATE rawhash_ruclient.cpp)
        add_dependencies(${TARGET_NAME} ruclient_build)
    endif()
endfunction()

function(setup_ruclient)
    if(RUCLIENT_ENABLED)
        if(NOT RUCLIENT_DIR)
            override_cached(RUCLIENT_DIR ${EXTERNAL_PROJECTS_BUILD_DIR}/ruclient)
        endif()
        ExternalProject_Add(
            ruclient_build
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/readuntil_fake
            BINARY_DIR ${RUCLIENT_DIR}/build
            CMAKE_ARGS
                -DCMAKE_INSTALL_PREFIX=${RUCLIENT_DIR}
        )
        include_directories(${RUCLIENT_DIR}/include)
        message(STATUS "ruclient enabled")
    else()
        message(STATUS "ruclient disabled")
    endif()
endfunction()
