include(ExternalProject)

function(override_cached name value)
    message(STATUS "Overriding ${name} to ${value}")
    get_property(doc_string CACHE ${name} PROPERTY HELPSTRING)
    get_property(type CACHE ${name} PROPERTY TYPE)
    set(${name} "${value}" CACHE ${type} "${doc_string}" FORCE)
endfunction()

function(link_imported_library LIB_NAME LIB_DIR)
    add_library(${LIB_NAME} STATIC IMPORTED)
    file(MAKE_DIRECTORY ${LIB_DIR}/include)
    set_target_properties(${LIB_NAME} PROPERTIES
        IMPORTED_LOCATION ${LIB_DIR}/lib/lib${LIB_NAME}.a
        INTERFACE_INCLUDE_DIRECTORIES ${LIB_DIR}/include
    )
    target_link_libraries(${TARGET_NAME} PRIVATE ${LIB_NAME})
endfunction()
