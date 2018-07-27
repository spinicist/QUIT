option( BUILD_RELAX "Build the relaxometry (DESPOT etc.) programs" ON )
if( ${BUILD_RELAX} )
    set( PROGRAMS
        qiafi qidream
        qidespot1 qidespot2 qidespot1hifi qidespot2fm
        qimultiecho qimp2rage
       )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_sequences qi_imageio qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    # mcdespot is a special flower
    add_executable( qimcdespot qimcdespot.cpp TwoPoolModel.cpp ThreePoolModel.cpp Helpers.cpp )
    target_link_libraries( qimcdespot qi_sequences qi_imageio qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES} )

    install( TARGETS ${PROGRAMS} qimcdespot RUNTIME DESTINATION bin )
endif()