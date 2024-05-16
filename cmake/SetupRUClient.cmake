function(setup_ruclient)
    if(RUCLIENT_ENABLED)
        find_package(gRPC REQUIRED)
        set_target_properties(${TARGET_NAME} PROPERTIES CXX_STANDARD 20)
        set_target_properties(${TARGET_NAME} PROPERTIES C_STANDARD 20)
        target_compile_definitions(${TARGET_NAME} PRIVATE RUCLIENT_ENABLED)
        target_sources(${TARGET_NAME} PRIVATE rawhash_ruclient.cpp)

        set(RUCLIENT_DIR ${CMAKE_SOURCE_DIR}/extern/readuntil_fake)
        
        add_subdirectory(${RUCLIENT_DIR} ${WORKDIR}/ruclient EXCLUDE_FROM_ALL)
        message(STATUS "ruclient enabled")
    else()
        message(STATUS "ruclient disabled")
    endif()
endfunction()
