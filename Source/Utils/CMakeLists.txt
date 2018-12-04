option( BUILD_UTILITIES "Build the utility programs" ON )
if( ${BUILD_UTILITIES} )
    set( PROGRAMS
            qiaffine
            qi_coil_combine
            qicomplex
            qi_gradient
            qikfilter    
            qimask
            qi_pca
            qipolyfit
            qipolyimg
            qi_rfprofile
            qi_tvmask
    )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_imageio qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
endif()