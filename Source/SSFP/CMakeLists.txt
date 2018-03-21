option( BUILD_ELLIPSE "Build the SSFP programs" ON )
if( ${BUILD_ELLIPSE} )
    add_executable( qi_ssfp_bands qi_ssfp_bands.cpp Banding.cpp )
    target_link_libraries( qi_ssfp_bands qi_sequences qi_imageio qi_filters qi_core )

    add_executable( qi_ssfp_ellipse qi_ssfp_ellipse.cpp
        EllipseAlgo.cpp EllipseHelpers.cpp
        DirectAlgo.cpp Direct2Algo.cpp HyperAlgo.cpp )
    target_link_libraries( qi_ssfp_ellipse qi_sequences qi_imageio qi_filters qi_core ${CERES_LIBRARIES} )

    add_executable( qi_ssfp_planet qi_ssfp_planet.cpp )
    target_link_libraries( qi_ssfp_planet qi_sequences qi_imageio qi_filters qi_core )

    add_executable( qi_ssfp_emt qi_ssfp_emt.cpp MTFromEllipse.cpp )
    target_link_libraries( qi_ssfp_emt qi_sequences qi_imageio qi_filters qi_core ${CERES_LIBRARIES} )

    install( TARGETS qi_ssfp_bands qi_ssfp_ellipse qi_ssfp_planet qi_ssfp_emt RUNTIME DESTINATION bin )
endif()