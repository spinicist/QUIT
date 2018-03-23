option( BUILD_CEST "Build the MT Programs" ON )
if( ${BUILD_CEST} )
    set( PROGRAMS
         qi_mtasym qi_lorentzian qi_dipolar_mtr )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_sequences qi_filters qi_imageio qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
endif()