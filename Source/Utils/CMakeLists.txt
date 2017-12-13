set( PROGRAMS
    qihdr qinewimage qidiff
    qicomplex qiaffine qimask qikfilter qireorder
    qisplitsubjects qiunwrap qipolyfit qipolyimg
    qi_coil_combine qi_bet qi_rfprofile
    qisignal )

foreach(PROGRAM ${PROGRAMS})
    add_executable(${PROGRAM} ${PROGRAM}.cpp)
    target_link_libraries(${PROGRAM} qi_sequences qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
endforeach(PROGRAM)

install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )