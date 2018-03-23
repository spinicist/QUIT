set( PROGRAMS qinewimage qisignal qidiff )

foreach( PROGRAM ${PROGRAMS} )
    add_executable( ${PROGRAM} ${PROGRAM}.cpp )
    target_link_libraries( ${PROGRAM} qi_sequences qi_filters qi_imageio qi_core ${ITK_LIBRARIES} )
endforeach( PROGRAM )

install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )