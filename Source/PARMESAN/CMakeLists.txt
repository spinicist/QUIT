option( BUILD_PARMESAN "Build the relaxometry (DESPOT etc.) programs" OFF )
if( ${BUILD_PARMESAN} )
    set(SOURCES
            prep_sequence.cpp
            parmesan_main.cpp
            GridInterp.cpp
            prep_qmt_Rx.cpp
            prep_qmt_k.cpp
        )
    target_sources(qi PRIVATE ${SOURCES})
    target_compile_definitions(qi PRIVATE "-DBUILD_PARMESAN")
endif()