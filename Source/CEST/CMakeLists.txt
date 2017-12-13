option( BUILD_CEST "Build the CEST Programs" ON )
if( ${BUILD_CEST} )
    set( PROGRAMS
        qi_cestasym qi_cestlorentz )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_sequences qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
endif()