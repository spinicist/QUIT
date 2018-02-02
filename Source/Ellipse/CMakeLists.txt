option( BUILD_ELLIPSE "Build the SSFP Ellipse programs" ON )
if( ${BUILD_ELLIPSE} )
    add_executable( qissfpbands qissfpbands.cpp Banding.cpp )
    target_link_libraries( qissfpbands qi_sequences qi_imageio qi_filters qi_core ${ITK_LIBRARIES} )

    add_executable( qi_ellipse qi_ellipse.cpp
        EllipseAlgo.cpp EllipseHelpers.cpp
        DirectAlgo.cpp Direct2Algo.cpp HyperAlgo.cpp )
    target_link_libraries( qi_ellipse qi_sequences qi_imageio qi_filters qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES} )

    add_executable( qi_emt qi_emt.cpp MTFromEllipse.cpp )
    target_link_libraries( qi_emt qi_sequences qi_imageio qi_filters qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES} )

    install( TARGETS qissfpbands qi_ellipse qi_emt RUNTIME DESTINATION bin )
endif()