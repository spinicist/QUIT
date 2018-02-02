option( BUILD_RELAX "Build the relaxometry (DESPOT etc.) programs" ON )
if( ${BUILD_RELAX} )
    set( PROGRAMS
        qiafi qidream
        qidespot1 qidespot2 qidespot1hifi qidespot2fm qimcdespot 
        qimultiecho qimp2rage )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_sequences qi_imageio qi_filters qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
endif()