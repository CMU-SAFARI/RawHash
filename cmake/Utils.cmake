include(ExternalProject)

function(override_cached name value)
    message(STATUS "Overriding ${name} to ${value}")
    get_property(doc_string CACHE ${name} PROPERTY HELPSTRING)
    get_property(type CACHE ${name} PROPERTY TYPE)
    set(${name} ${value} CACHE ${type} ${doc_string} FORCE)
endfunction()


function(add_imported_library TARGET_NAME LIB_NAME)
    target_link_libraries(${TARGET_NAME} PRIVATE ${LIB_NAME}_RAWHASH_LIB)
endfunction()

function(define_imported_library LIB_NAME LIB_DIR LIB_TYPE)
    add_library(${LIB_NAME}_RAWHASH_LIB ${LIB_TYPE} IMPORTED)
    set(EXTENSION ".a")
    if(LIB_TYPE STREQUAL SHARED)
        set(EXTENSION ".so")
    endif()
    set_target_properties(${LIB_NAME}_RAWHASH_LIB PROPERTIES
        IMPORTED_LOCATION ${LIB_DIR}/lib/lib${LIB_NAME}${EXTENSION}
        INTERFACE_INCLUDE_DIRECTORIES ${LIB_DIR}/include)
    file(MAKE_DIRECTORY ${LIB_DIR}/include)
    # Can't install(TARGETS ...) for external projects
    # Also some .so are symlinks, so install all
    install(DIRECTORY ${LIB_DIR}/lib/ DESTINATION lib
            FILES_MATCHING PATTERN *${EXTENSION}*)
    install(DIRECTORY ${LIB_DIR}/include/
            DESTINATION include/${PROJECT_NAME})
endfunction()


function(check_directory_exists_and_non_empty DIR)
    if(NOT EXISTS ${DIR})
        message(FATAL_ERROR "Directory ${DIR} does not exist.")
    endif()

    file(GLOB DIRECTORY_CONTENTS "${DIR}/*")
    if(DIRECTORY_CONTENTS STREQUAL "")
        message(FATAL_ERROR "Directory ${DIR} is empty.")
    endif()

    message(STATUS "Directory ${DIR} exists and is non-empty.")
endfunction()
