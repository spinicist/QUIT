option( BUILD_PERFUSION "Build the perfusion programs" ON )
if( ${BUILD_PERFUSION} )
    set( PROGRAMS qi_asl qi_ase_oef qi_zshim )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_imageio qi_sequences qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
endif()