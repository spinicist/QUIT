option( BUILD_PARMESAN "Build the relaxometry (DESPOT etc.) programs" OFF )
if( ${BUILD_PARMESAN} )
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/prep_sequence.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/parmesan_main.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/GridInterp.cpp
        )
    target_sources(qi PRIVATE ${SOURCES})
    target_compile_definitions(qi PRIVATE "-DBUILD_PARMESAN")
endif()