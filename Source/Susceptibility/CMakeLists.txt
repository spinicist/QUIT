option( BUILD_SUSCEPTIBILITY "Build the susceptibility programs" ON )
if( ${BUILD_SUSCEPTIBILITY} )
    add_executable(qi_unwrap_laplace qi_unwrap_laplace.cpp)
    target_link_libraries(qi_unwrap_laplace qi_imageio qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})

    add_executable( qi_unwrap_path qi_unwrap_path.cpp
                        ReliabilityFilter.cpp
                        PathUnwrapFilter.cpp )
    target_link_libraries( qi_unwrap_path qi_imageio qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES} )

    add_executable(qi_fieldmap qi_fieldmap.cpp)
    target_link_libraries(qi_fieldmap qi_imageio qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})

    install( TARGETS qi_unwrap_laplace qi_unwrap_path qi_fieldmap RUNTIME DESTINATION bin )
endif()