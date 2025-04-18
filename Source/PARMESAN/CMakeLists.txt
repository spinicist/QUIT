option( BUILD_PARMESAN "Build the relaxometry (DESPOT etc.) programs" OFF )
if( ${BUILD_PARMESAN} )
    set(SOURCES
            Basis.cpp
            parmesan_main.cpp
            parmesan_sc.cpp
            GridInterp.cpp
            ParameterGrid.cpp
            prep_qmt_Rx.cpp
            prep_qmt_k.cpp
            prep_sc.cpp
            prep_sequence.cpp
        )
    target_sources(qi PRIVATE ${SOURCES})
    target_compile_definitions(qi PRIVATE "-DBUILD_PARMESAN")
endif()