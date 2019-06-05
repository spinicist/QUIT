option( BUILD_MT "Build the MT Programs" ON )
if( ${BUILD_MT} )
    set( PROGRAMS
         qi_zspec_interp qi_zspec_b1 qi_lorentzian
         qi_lineshape qi_qmt qi_mtsat
       )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_sequences qi_imageio qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
endif()