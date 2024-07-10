include(${CMAKE_CURRENT_LIST_DIR}/Util.cmake)

function(add_tflite_to_target TARGET_NAME)
    target_link_libraries(${TARGET_NAME} PRIVATE tensorflow-lite)
endfunction()

function(setup_tflite)
    set(TF_SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/tensorflow)
    add_subdirectory(${TF_SOURCE_DIR}/tensorflow/lite ${EXTERNAL_PROJECTS_BUILD_DIR}/tflite EXCLUDE_FROM_ALL)
    include_directories(${TF_SOURCE_DIR})
endfunction()
