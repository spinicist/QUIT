option( BUILD_RUFIS "Build the relaxometry (DESPOT etc.) programs" ON )
if( ${BUILD_RUFIS} )
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/rf_pulse.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/rf_sim.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ss_sequence.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ss_model.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ss_T2.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ss_mt.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ss_main.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/transient_sequence.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/transient_b1_model.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/transient_mt_model.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/transient_main.cpp
        )
    target_sources(qi PRIVATE ${SOURCES})
    target_compile_definitions(qi PRIVATE "-DBUILD_RUFIS")
endif()