option( BUILD_PARMESAN "Build the relaxometry (DESPOT etc.) programs" OFF )
if( ${BUILD_PARMESAN} )
    set(SOURCES
            Basis.cpp
            GridInterp.cpp
            ParameterGrid.cpp
            parmesan_sc.cpp
            parmesan_opt.cpp
            parmesan_qmt.cpp
            prep_sc.cpp
            prep_qmt.cpp
            prep_sequence.cpp
        )
    target_sources(qi PRIVATE ${SOURCES})
    target_compile_definitions(qi PRIVATE "-DBUILD_PARMESAN")
    find_package(pagmo REQUIRED)
    target_link_libraries(qi PRIVATE
        Pagmo::pagmo
    )
endif()