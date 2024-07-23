include(${CMAKE_CURRENT_LIST_DIR}/Utils.cmake)

function(add_zstd_to_target TARGET_NAME)
    set(ZSTD_DIR ${EXTERNAL_PROJECTS_BUILD_DIR}/zstd)
    add_dependencies(${TARGET_NAME} zstd_build)
    target_link_libraries(${TARGET_NAME} PRIVATE zstd)
endfunction()

function(setup_zstd)
    set(ZSTD_DIR ${EXTERNAL_PROJECTS_BUILD_DIR}/zstd)
    ExternalProject_Add(
        zstd_build
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/zstd/build/cmake
        BINARY_DIR ${ZSTD_DIR}/build
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${ZSTD_DIR}
    )
    define_imported_library(zstd ${ZSTD_DIR})
endfunction()

function(add_pod5_to_target TARGET_NAME)
    if(NOPOD5)
        target_compile_definitions(${TARGET_NAME} PRIVATE NPOD5RH=1)
    else()
        add_zstd_to_target(${TAARGET_NAME})

        set(POD5_VERSION "0.2.2")
        set(POD5_URLDIR "pod5-${POD5_VERSION}-${CMAKE_SYSTEM_NAME}")
        set(POD5_REPO "https://github.com/nanoporetech/pod5-file-format")

        resolve_pod5_url()

        if(POD5_DOWNLOAD)
            add_dependencies(${TARGET_NAME} pod5_download)
        endif()
        target_link_libraries(${TARGET_NAME} PRIVATE ${POD5_LIBRARIES} zstd)
    endif()
endfunction()

function(setup_pod5)
    if(NOT NOPOD5)
        setup_zstd()

        set(POD5_VERSION "0.2.2")
        set(POD5_URLDIR "pod5-${POD5_VERSION}-${CMAKE_SYSTEM_NAME}")
        set(POD5_REPO "https://github.com/nanoporetech/pod5-file-format")

        resolve_pod5_url()

        if(POD5_DOWNLOAD)
            if(NOT POD5_DIR)
                override_cached(POD5_DIR ${EXTERNAL_PROJECTS_BUILD_DIR}/${POD5_URLDIR})
            endif()
            ExternalProject_Add(
                pod5_download
                SOURCE_DIR ${POD5_DIR}
                URL ${POD5_URL}
                CONFIGURE_COMMAND ""
                BUILD_COMMAND ""
                INSTALL_COMMAND ""
                # requires cmake >= 3.24
                # DOWNLOAD_EXTRACT_TIMESTAMP TRUE
            )
            add_dependencies(${TARGET_NAME} pod5_download)
        else()
            if(NOT POD5_DIR)
                message(FATAL_ERROR "POD5_DOWNLOAD is OFF, but no dir provided")
            endif()
        endif()
        include_directories(${POD5_DIR}/include)
    endif()
endfunction()

# POD5_URL and POD5_LIBRARIES are set at PARENT_SCOPE
function(resolve_pod5_url)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(POD5_LIB "lib64")
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm)")
            if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
                set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-gcc7-arm64.tar.gz")
                set(POD5_LIB "lib")
            else()
                set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-arm64.tar.gz" PARENT_SCOPE)
            endif()
        else()
            set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-x64.tar.gz" PARENT_SCOPE)
        endif()
        set(POD5_LIB_DIR "${EXTERNAL_PROJECTS_BUILD_DIR}/${POD5_URLDIR}/${POD5_LIB}")
        set(POD5_LIBRARIES "${POD5_LIB_DIR}/libpod5_format.so"
                        "${POD5_LIB_DIR}/libarrow.so"
                        "${POD5_LIB_DIR}/libjemalloc_pic.so"  PARENT_SCOPE)
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-osx-11.0-arm64.tar.gz")
        set(POD5_LIB_DIR "${EXTERNAL_PROJECTS_BUILD_DIR}/${POD5_URLDIR}/lib")
        set(POD5_LIBRARIES "${POD5_LIB_DIR}/libpod5_format.so"
                           "${POD5_LIB_DIR}/libarrow.so"  PARENT_SCOPE)
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows_NT")
        set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-win" PARENT_SCOPE)
    endif()
endfunction()


# not working because of improper design, PARENT_SCOPE should not be used, rather define targets properly
# include(${CMAKE_CURRENT_LIST_DIR}/Utils.cmake)

# set(ZSTD_DIR ${EXTERNAL_PROJECTS_BUILD_DIR}/zstd)

# function(setup_zstd)
#     ExternalProject_Add(
#         zstd_build
#         SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/zstd/build/cmake
#         BINARY_DIR ${ZSTD_DIR}/build
#         CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${ZSTD_DIR}
#     )
# endfunction()
    
# function(add_zstd_to_target TARGET_NAME)
#     add_dependencies(${TARGET_NAME} zstd_build)
#     link_imported_library(${TARGET_NAME} zstd ${ZSTD_DIR})
#     target_link_libraries(${TARGET_NAME} PRIVATE zstd)
# endfunction()

# function(setup_pod5)
#     if(NOPOD5)
#         return()
#     endif()

#     setup_zstd()

#     set(POD5_VERSION "0.2.2")
#     set(POD5_URLDIR "pod5-${POD5_VERSION}-${CMAKE_SYSTEM_NAME}")
#     set(POD5_REPO "https://github.com/nanoporetech/pod5-file-format")

#     resolve_pod5_url()

#     if(POD5_DOWNLOAD)
#         if(NOT POD5_DIR)
#             override_cached(POD5_DIR ${EXTERNAL_PROJECTS_BUILD_DIR}/${POD5_URLDIR})
#         endif()
#         ExternalProject_Add(
#             pod5_download
#             SOURCE_DIR ${POD5_DIR}
#             URL ${POD5_URL}
#             CONFIGURE_COMMAND ""
#             BUILD_COMMAND ""
#             INSTALL_COMMAND ""
#             # requires cmake >= 3.24
#             # DOWNLOAD_EXTRACT_TIMESTAMP TRUE
#         )
#     else()
#         if(NOT POD5_DIR)
#             message(FATAL_ERROR "POD5_DOWNLOAD is OFF, but no dir provided")
#         endif()
#     endif()
# endfunction()

# function(add_pod5_to_target TARGET_NAME)
#     if(NOPOD5)
#         target_compile_definitions(${TARGET_NAME} PRIVATE NPOD5RH=1)
#     else()
#         add_zstd_to_target(${TARGET_NAME})

#         add_dependencies(${TARGET_NAME} pod5_download)
#         include_directories(${POD5_DIR}/include)
#         message(STATUS "Adding include dir ${POD5_DIR}/include, POD5 libraries: ${POD5_LIBRARIES}")
#         target_link_libraries(${TARGET_NAME} PRIVATE ${POD5_LIBRARIES})
#     endif()
# endfunction()

# # POD5_URL and POD5_LIBRARIES are set at PARENT_SCOPE
# function(resolve_pod5_url)
#     if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
#         set(POD5_LIB "lib64")
#         if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm)")
#             if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
#                 set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-gcc7-arm64.tar.gz")
#                 set(POD5_LIB "lib")
#             else()
#                 set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-arm64.tar.gz" PARENT_SCOPE)
#             endif()
#         else()
#             set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-x64.tar.gz" PARENT_SCOPE)
#         endif()
#         set(POD5_LIB_DIR "${EXTERNAL_PROJECTS_BUILD_DIR}/${POD5_URLDIR}/${POD5_LIB}")
#         set(POD5_LIBRARIES "${POD5_LIB_DIR}/libpod5_format.a"
#                         "${POD5_LIB_DIR}/libarrow.a"
#                         "${POD5_LIB_DIR}/libjemalloc_pic.a"  PARENT_SCOPE)
#     elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
#         set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-osx-11.0-arm64.tar.gz")
#         set(POD5_LIB_DIR "${EXTERNAL_PROJECTS_BUILD_DIR}/${POD5_URLDIR}/lib")
#         set(POD5_LIBRARIES "${POD5_LIB_DIR}/libpod5_format.a"
#                            "${POD5_LIB_DIR}/libarrow.a"  PARENT_SCOPE)
#     elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows_NT")
#         set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-win" PARENT_SCOPE)
#         # todo: not setting libraries!
#     endif()
# endfunction()