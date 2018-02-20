option( BUILD_SUSCEPTIBILITY "Build the susceptibility programs" ON )
if( ${BUILD_SUSCEPTIBILITY} )
    set( PROGRAMS qi_unwrap_laplace qi_unwrap_path )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_imageio qi_filters qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
endif()