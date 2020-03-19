option( BUILD_RUFIS "Build the relaxometry (DESPOT etc.) programs" ON )
if( ${BUILD_RUFIS} )
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/Commands.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/rufis_pulse.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/rufis_sequence.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/rf_sim.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ss_sequence.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ss_main.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ss_model.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ss_T2.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/mupa_main.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/mupa_model_b1.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/mupa_model_mt.cpp
        )
    target_sources(qi PRIVATE ${SOURCES})
    target_compile_definitions(qi PRIVATE "-DBUILD_RUFIS")
endif()