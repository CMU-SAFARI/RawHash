include(${CMAKE_CURRENT_LIST_DIR}/Util.cmake)

function(setup_hdf5 TARGET_NAME)
    if(NOHDF5)
        target_compile_definitions(${TARGET_NAME} PRIVATE NHDF5RH=1)
    else()
        if(HDF5_COMPILE)
            if(NOT HDF5_DIR)
                override_cached(HDF5_DIR ${WORKDIR}/hdf5)
            endif()
            set(HDF5_BUILD_DIR ${HDF5_DIR}/build)
            ExternalProject_Add(
                hdf5_build
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/hdf5
                BINARY_DIR ${HDF5_BUILD_DIR}
                CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/extern/hdf5/configure --enable-threadsafe --disable-hl --prefix=${HDF5_BUILD_DIR}
                # INSTALL_DIR and DCMAKE_INSTALL_PREFIX are ignored by hdf5
                INSTALL_COMMAND make install prefix=${HDF5_DIR}
            )
            add_dependencies(${TARGET_NAME} hdf5_build)
        else()
            if(NOT HDF5_DIR)
                message(FATAL_ERROR "HDF5_COMPILE is OFF, but no dir provided")
            endif()
        endif()
        link_imported_library(${TARGET_NAME} hdf5 ${HDF5_DIR})
    endif()
endfunction()
