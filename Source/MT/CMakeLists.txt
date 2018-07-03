option( BUILD_MT "Build the MT Programs" ON )
if( ${BUILD_MT} )
    set( PROGRAMS
        #  qi_mtasym qi_lorentzian qi_dipolar_mtr
         qi_lineshape )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
endif()