option( BUILD_RELAX "Build the relaxometry (DESPOT etc.) programs" ON )
if( ${BUILD_RELAX} )
    set( PROGRAMS
        qidespot1 qidespot2 qidespot1hifi qidespot2fm qi_jsr
        qimultiecho qimp2rage qi_mpm_r2s
       )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_sequences qi_imageio qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    # mcdespot is a special flower
    # add_executable( qimcdespot qimcdespot.cpp TwoPoolModel.cpp ThreePoolModel.cpp Helpers.cpp )
    # target_link_libraries( qimcdespot qi_sequences qi_imageio qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES} )

    install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
endif()