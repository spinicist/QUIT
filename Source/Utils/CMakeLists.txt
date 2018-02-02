option( BUILD_UTILITIES "Build the utility programs" ON )
if( ${BUILD_UTILITIES} )
    set( PROGRAMS
        qihdr
        qicomplex qiaffine qimask qikfilter qireorder
        qisplitsubjects qiunwrap qipolyfit qipolyimg
        qi_coil_combine qi_bet qi_rfprofile )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_imageio qi_filters qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
endif()