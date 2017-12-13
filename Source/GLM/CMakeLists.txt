set( PROGRAMS
    qi_rois qi_glmsetup qi_glmcontrasts )

foreach(PROGRAM ${PROGRAMS})
    add_executable(${PROGRAM} ${PROGRAM}.cpp)
    target_link_libraries(${PROGRAM} qi_sequences qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
endforeach(PROGRAM)

install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )