include(${CMAKE_CURRENT_LIST_DIR}/Util.cmake)

function(add_hdf5_to_target TARGET_NAME)
    if(NOHDF5)
        target_compile_definitions(${TARGET_NAME} PRIVATE NHDF5RH=1)
    else()
        if(HDF5_COMPILE)
            add_dependencies(${TARGET_NAME} hdf5_build)
        endif()
        link_imported_library(${TARGET_NAME} hdf5 ${HDF5_DIR})
    endif()
endfunction()

function(setup_hdf5)
    if(NOT NOHDF5)
        # print HDF5_DIR
        message(STATUS "EXTERNAL_PROJECTS_BUILD_DIR: ${EXTERNAL_PROJECTS_BUILD_DIR}")
        message(STATUS "HDF5_DIR: ${HDF5_DIR}")
        set(HDF5_SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/hdf5)
        if(HDF5_COMPILE)
            if(NOT HDF5_DIR)
                override_cached(HDF5_DIR ${EXTERNAL_PROJECTS_BUILD_DIR}/hdf5)
            endif()
            set(HDF5_BUILD_DIR ${HDF5_DIR}/build)
            ExternalProject_Add(
                hdf5_build
                SOURCE_DIR ${HDF5_SOURCE_DIR}
                BINARY_DIR ${HDF5_BUILD_DIR}
                CONFIGURE_COMMAND ${HDF5_SOURCE_DIR}/configure --enable-threadsafe --disable-hl --prefix=${HDF5_BUILD_DIR}
                # INSTALL_DIR and DCMAKE_INSTALL_PREFIX are ignored by hdf5
                INSTALL_COMMAND make install prefix=${HDF5_DIR}
            )
        else()
            if(NOT HDF5_DIR)
                message(FATAL_ERROR "HDF5_COMPILE is OFF, but no dir provided")
            endif()
        endif()
        define_imported_library(hdf5 ${HDF5_DIR})
    endif()
endfunction()
