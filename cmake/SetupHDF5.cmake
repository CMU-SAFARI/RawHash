include(ExternalProject)

function(setup_hdf5)
    if(NOHDF5)
        target_compile_definitions(${TARGET_NAME} PRIVATE NHDF5RH=1)
    else()
        set(HDF5_DIR "${WORKDIR}/hdf5")
        set(HDF5_BUILD_DIR "${HDF5_DIR}/build")
        if(HDF5_COMPILE)
            if(NOT HDF5_INCLUDE_DIR)
                override_cached(HDF5_INCLUDE_DIR "${HDF5_DIR}/include")
            endif()
            ExternalProject_Add(
                hdf5_build
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/hdf5
                BINARY_DIR ${HDF5_BUILD_DIR}
                CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/extern/hdf5/configure --enable-threadsafe --disable-hl --prefix=${HDF5_BUILD_DIR}
                INSTALL_COMMAND make install prefix=${HDF5_DIR}
            )
            add_dependencies(${TARGET_NAME} hdf5_build)
        else()
            if(NOT HDF5_INCLUDE_DIR)
                message(FATAL_ERROR "HDF5_COMPILE is OFF, but no include dir provided")
            endif()
        endif()
        add_library(hdf5 STATIC IMPORTED)
        file(MAKE_DIRECTORY ${HDF5_DIR}/include)
        set_target_properties(hdf5 PROPERTIES
            IMPORTED_LOCATION $${HDF5_DIR}/lib/libhdf5.a
            INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIR}
        )
        target_link_libraries(${TARGET_NAME} PRIVATE hdf5)
    endif()
endfunction()
