option( BUILD_PARMESAN "Build the relaxometry (DESPOT etc.) programs" OFF )
if( ${BUILD_PARMESAN} )
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/rf_pulse.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/rf_sim.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/transient_sequence.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/transient_qmt_model.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/transient_main.cpp
        )
    target_sources(qi PRIVATE ${SOURCES})
    target_compile_definitions(qi PRIVATE "-DBUILD_PARMESAN")
endif()