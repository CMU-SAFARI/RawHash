include(ExternalProject)

function(override_cached name value)
    message(STATUS "Overriding ${name} to ${value}")
    get_property(doc_string CACHE ${name} PROPERTY HELPSTRING)
    get_property(type CACHE ${name} PROPERTY TYPE)
    set(${name} ${value} CACHE ${type} ${doc_string} FORCE)
endfunction()

function(link_imported_library TARGET_NAME LIB_NAME LIB_DIR)
    add_library(${LIB_NAME} SHARED IMPORTED)
    file(MAKE_DIRECTORY ${LIB_DIR}/include)
    set_target_properties(${LIB_NAME} PROPERTIES
        IMPORTED_LOCATION ${LIB_DIR}/lib/lib${LIB_NAME}.so
        INTERFACE_INCLUDE_DIRECTORIES ${LIB_DIR}/include
    )
    target_link_libraries(${TARGET_NAME} PRIVATE ${LIB_NAME})
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
