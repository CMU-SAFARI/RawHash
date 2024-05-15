include(${CMAKE_CURRENT_LIST_DIR}/Util.cmake)

function(setup_tflite)
    set(TF_SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/tensorflow)
    add_subdirectory(${TF_SOURCE_DIR}/tensorflow/lite ${WORKDIR}/tflite EXCLUDE_FROM_ALL)
    include_directories(${TF_SOURCE_DIR})
    target_link_libraries(${TARGET_NAME} PRIVATE tensorflow-lite)
endfunction()
